#!/usr/bin/env python3
"""
gel_autoalign.py - automatically straighten and crop gel images (TIFF).

What it does
------------
1. Finds the gel slab (the gel-vs-background contrast) with classical image
   processing: Gaussian blur -> Otsu threshold -> largest object -> contour.
2. Measures the tilt by fitting straight lines to the four edges of the slab
   (initialised from OpenCV's minimum-area rectangle).
3. Rotates by exactly that angle and crops to the slab plus a small pad, in a
   SINGLE resampling step. If the tilt is negligible, it does a pure integer
   crop with no interpolation at all.
4. Writes a lossless TIFF (same bit depth, channels and DPI as the input) and
   a JSON sidecar recording everything needed to reproduce the result.

Deterministic: no machine learning, no random numbers. Same input + same
parameters + same library versions -> byte-identical output image.
Originals are never modified; outputs go to a separate folder.

Usage
-----
    python gel_autoalign.py gel.tif                  # one file
    python gel_autoalign.py /path/to/folder          # every .tif/.tiff in it
    python gel_autoalign.py /path/to/folder --watch  # keep processing new files
    python gel_autoalign.py gel.tif --pad 40 --debug # bigger margin + QC image
    python gel_autoalign.py gel.tif --label "MP-02 PCR"   # also write a labeled copy
    python gel_autoalign.py folder --label "{name}"       # label = each file's name

Environment (conda/mamba):
    mamba env create -f environment.yml
    mamba activate gel-autoalign

Labeling: --label writes an ADDITIONAL 8-bit RGB display copy (<name>_aligned_labeled.tif)
with the text burned in at the gel's bottom-right, set in Latin Modern Sans Bold
(bundled in fonts/). The plain <name>_aligned.tif is never
modified by labeling, so it stays the clean record.

Orientation note: this removes tilt (modulo 90 degrees). It does not decide
which side is "up", so keep placing gels on the imager in the same orientation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
import tomllib
from datetime import datetime, timezone
from pathlib import Path

import cv2
import numpy as np
import tifffile

__version__ = "1.2.2"

INTERP = {
    "cubic": cv2.INTER_CUBIC,
    "linear": cv2.INTER_LINEAR,
    "lanczos": cv2.INTER_LANCZOS4,
}


# --------------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tiff(path: Path):
    """Return (array HxW or HxWxC, resolution dict or None, notes list)."""
    notes = []
    with tifffile.TiffFile(path) as tf:
        if len(tf.pages) > 1:
            notes.append(f"file has {len(tf.pages)} pages; only the first was used")
        page = tf.pages[0]
        arr = page.asarray()
        res = None
        try:
            xr = page.tags.get("XResolution")
            yr = page.tags.get("YResolution")
            ru = page.tags.get("ResolutionUnit")
            if xr is not None and yr is not None:
                n, d = xr.value
                m, e = yr.value
                res = {
                    "resolution": (n / d, m / e),
                    "resolutionunit": int(ru.value) if ru is not None else 2,
                }
        except Exception:
            res = None

    arr = np.squeeze(arr)
    if arr.ndim == 3 and arr.shape[0] in (3, 4) and arr.shape[-1] not in (3, 4):
        arr = np.moveaxis(arr, 0, -1)  # planar -> interleaved
    if arr.ndim == 3 and arr.shape[-1] == 4:
        notes.append("alpha channel dropped")
        arr = arr[..., :3]
    if arr.ndim not in (2, 3) or (arr.ndim == 3 and arr.shape[-1] != 3):
        raise ValueError(f"unsupported image shape {arr.shape}")
    if arr.dtype not in (np.uint8, np.uint16, np.float32):
        raise ValueError(f"unsupported dtype {arr.dtype} (use 8/16-bit or float32)")
    return np.ascontiguousarray(arr), res, notes


def git_info():
    """Git state of this script's own checkout, or None if not in a git checkout."""
    import subprocess

    repo = Path(__file__).resolve().parent

    def git(*cmd):
        return subprocess.run(
            ["git", "-C", str(repo), *cmd], capture_output=True, text=True, check=True, timeout=10
        ).stdout.strip()

    try:
        commit = git("rev-parse", "HEAD")
        describe = git("describe", "--tags", "--always", "--dirty")
        try:
            remote = git("config", "--get", "remote.origin.url") or None
        except subprocess.CalledProcessError:
            remote = None
        return {
            "describe": describe,
            "commit": commit,
            "remote": remote,
            "dirty": describe.endswith("-dirty"),
        }
    except Exception:
        return None


def write_tiff(path: Path, img: np.ndarray, res, description: str):
    kw = {}
    if res:
        kw["resolution"] = res["resolution"]
        kw["resolutionunit"] = res["resolutionunit"]
    photometric = "rgb" if img.ndim == 3 else "minisblack"
    tifffile.imwrite(
        path,
        img,
        photometric=photometric,
        compression="zlib",  # lossless
        description=description,
        metadata=None,
        **kw,
    )


# --------------------------------------------------------------------------
# Gel detection
# --------------------------------------------------------------------------
def to_gray(arr: np.ndarray) -> np.ndarray:
    a = arr.astype(np.float32)
    return a.mean(axis=2) if a.ndim == 3 else a


def border_background(arr: np.ndarray):
    """Median value of the outer image border, per channel (assumes the border is background)."""
    h, w = arr.shape[:2]
    b = max(2, int(0.01 * min(h, w)))
    parts = [arr[:b], arr[-b:], arr[:, :b], arr[:, -b:]]
    if arr.ndim == 3:
        flat = np.concatenate([p.reshape(-1, arr.shape[2]) for p in parts])
        return np.median(flat, axis=0)
    flat = np.concatenate([p.ravel() for p in parts])
    return np.array([np.median(flat)])


def _wrap90(a):
    """Wrap an angle in degrees into [-45, 45)."""
    return (a + 45.0) % 90.0 - 45.0


def _wrap180(a):
    """Wrap an angle in degrees into [-90, 90)."""
    return (a + 90.0) % 180.0 - 90.0


def detect_gel(arr: np.ndarray):
    """
    Locate the gel slab and return a dict with:
      phi        tilt of the slab in degrees (image coords, y down; + = clockwise)
      center     (cx, cy) of the slab in pixel coordinates
      size       (width, height) of the slab after de-rotation
      corners    4x2 array of slab corners in the input image
      warnings   list of strings
    """
    warnings = []
    gray = to_gray(arr)
    h, w = gray.shape
    sigma = max(1.5, 0.002 * max(h, w))
    g = cv2.GaussianBlur(gray, (0, 0), sigma)

    # Background level from the image border; gel = anything that differs from it.
    # Works for both white-background (inverted) and dark-background images.
    b = max(2, int(0.01 * min(h, w)))
    border = np.concatenate([g[:b].ravel(), g[-b:].ravel(), g[:, :b].ravel(), g[:, -b:].ravel()])
    bg = float(np.median(border))
    diff = np.abs(g - bg)

    scale = float(np.percentile(diff, 99.5))
    if scale < 1e-6:
        raise ValueError("image is uniform; no gel found")
    d8 = np.clip(diff / scale * 255.0, 0, 255).astype(np.uint8)
    _, mask = cv2.threshold(d8, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    # Remove specks / dust, then keep the largest object.
    k = max(3, int(0.005 * max(h, w)) | 1)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k, k))
    mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN, kernel)
    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        raise ValueError("no gel-like object found")
    cnt = max(contours, key=cv2.contourArea)
    solid = np.zeros_like(mask)
    cv2.drawContours(solid, [cnt], -1, 255, thickness=cv2.FILLED)  # fills inner holes (bands, etc.)
    area = float(cv2.countNonZero(solid))

    pts = cnt.reshape(-1, 2).astype(np.float64)
    hull = cv2.convexHull(cnt).reshape(-1, 2).astype(np.float32)

    # Initial tilt from the minimum-area rectangle.
    rect = cv2.minAreaRect(hull)
    box = cv2.boxPoints(rect)
    e = box[1] - box[0]
    phi0 = _wrap90(math.degrees(math.atan2(e[1], e[0])))

    # Refine: robust line fits to each of the four sides, using the full contour.
    def frame(phi):
        c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
        center = pts.mean(axis=0)
        rel = pts - center
        u = rel[:, 0] * c + rel[:, 1] * s
        v = -rel[:, 0] * s + rel[:, 1] * c
        return u, v

    u, v = frame(phi0)
    umin, umax, vmin, vmax = u.min(), u.max(), v.min(), v.max()
    W, H = umax - umin, vmax - vmin
    band = 0.04
    trim = 0.10
    sides = {
        "top": ((v < vmin + band * H) & (u > umin + trim * W) & (u < umax - trim * W), 0.0),
        "bottom": ((v > vmax - band * H) & (u > umin + trim * W) & (u < umax - trim * W), 0.0),
        "left": ((u < umin + band * W) & (v > vmin + trim * H) & (v < vmax - trim * H), 90.0),
        "right": ((u > umax - band * W) & (v > vmin + trim * H) & (v < vmax - trim * H), 90.0),
    }
    devs, weights, side_devs = [], [], {}
    for name, (sel, base) in sides.items():
        p = pts[sel].astype(np.float32)
        if len(p) < 30:
            continue
        vx, vy, _, _ = cv2.fitLine(p, cv2.DIST_HUBER, 0, 0.01, 0.01).ravel()
        ang = math.degrees(math.atan2(vy, vx))
        dev = _wrap180(ang - (phi0 + base))
        devs.append(dev)
        weights.append(len(p))
        side_devs[name] = round(dev, 3)

    phi = phi0
    if devs:
        spread = max(devs) - min(devs)
        mean_dev = float(np.average(devs, weights=weights))
        if abs(mean_dev) > 2.0 or spread > 1.5:
            warnings.append(
                f"edge fits disagree (spread {spread:.2f} deg); used min-area rectangle angle instead"
            )
        else:
            phi = phi0 + mean_dev
    else:
        warnings.append("could not fit edges; used min-area rectangle angle")

    # Final tight extents at the refined angle.
    u, v = frame(phi)
    center0 = pts.mean(axis=0)
    c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    um, uM, vm, vM = u.min(), u.max(), v.min(), v.max()
    cu, cv_ = (um + uM) / 2, (vm + vM) / 2
    center = center0 + np.array([cu * c - cv_ * s, cu * s + cv_ * c])
    width, height = float(uM - um), float(vM - vm)
    hw, hh = width / 2, height / 2
    corners = np.array(
        [
            center + np.array([-hw * c + hh * s, -hw * s - hh * c]),
            center + np.array([hw * c + hh * s, hw * s - hh * c]),
            center + np.array([hw * c - hh * s, hw * s + hh * c]),
            center + np.array([-hw * c - hh * s, -hw * s + hh * c]),
        ]
    )

    # Quality checks.
    fill = area / max(width * height, 1.0)
    frac = area / (h * w)
    x, y, bw, bh = cv2.boundingRect(cnt)
    # Well notches (PAGE gels) legitimately remove ~10% of the outline near an edge,
    # (Mini-PROTEAN wells reach ~14% into the gel), so judge segmentation on the central
    # 60% x 60% of the slab, in the slab's own frame.
    R = np.array([[c, s], [-s, c]])
    cw, ch = max(int(math.ceil(width)), 1), max(int(math.ceil(height)), 1)
    t = np.array([(cw - 1) / 2, (ch - 1) / 2]) - R @ center
    core = cv2.warpAffine(solid, np.hstack([R, t.reshape(2, 1)]), (cw, ch), flags=cv2.INTER_NEAREST)
    core = core[int(0.2 * ch):int(0.8 * ch), int(0.2 * cw):int(0.8 * cw)]
    core_fill = float(np.count_nonzero(core)) / max(core.size, 1)
    if core_fill < 0.97:
        warnings.append(f"slab interior is only {core_fill:.0%} filled; segmentation may be incomplete")
    elif fill < 0.80:
        warnings.append(f"slab only fills {fill:.0%} of its bounding rectangle; segmentation may be incomplete")
    if frac < 0.02 or frac > 0.98:
        warnings.append(f"detected object covers {frac:.1%} of the image; check it is the gel")
    if x <= 1 or y <= 1 or x + bw >= w - 1 or y + bh >= h - 1:
        warnings.append("gel touches the image border; crop may be truncated")

    return {
        "phi": float(phi),
        "phi_minarearect": float(phi0),
        "side_deviations_deg": side_devs,
        "center": (float(center[0]), float(center[1])),
        "size": (width, height),
        "corners": corners,
        "fill_ratio": float(fill),
        "core_fill_ratio": core_fill,
        "area_fraction": float(frac),
        "background_gray": bg,
        "otsu_sigma_px": float(sigma),
        "warnings": warnings,
    }


# --------------------------------------------------------------------------
# Rotate + crop (single resampling step)
# --------------------------------------------------------------------------
def rotate_crop(arr, det, pad, interp, angle_tol):
    """Return (output image, info dict)."""
    phi = det["phi"]
    cx, cy = det["center"]
    width, height = det["size"]
    out_w = int(math.ceil(width + 2 * pad))
    out_h = int(math.ceil(height + 2 * pad))

    bg = border_background(arr)
    fill = tuple(float(x) for x in bg)
    ch = 1 if arr.ndim == 2 else arr.shape[2]
    fill = fill[:ch] if ch > 1 else (fill[0],)

    if abs(phi) < angle_tol:
        # Pure integer crop: exact pixel values, no interpolation.
        x0 = int(round(cx - (out_w - 1) / 2))
        y0 = int(round(cy - (out_h - 1) / 2))
        M = np.array([[1, 0, -x0], [0, 1, -y0]], dtype=np.float64)
        flags = cv2.INTER_NEAREST
        mode = "integer_crop_no_rotation"
    else:
        a = math.radians(phi)
        ca, sa = math.cos(a), math.sin(a)
        R = np.array([[ca, sa], [-sa, ca]], dtype=np.float64)  # rotates image by +phi CCW (visual)
        t = np.array([(out_w - 1) / 2, (out_h - 1) / 2]) - R @ np.array([cx, cy])
        M = np.hstack([R, t.reshape(2, 1)])
        flags = interp
        mode = "rotate_and_crop_single_resample"

    out = cv2.warpAffine(
        arr, M, (out_w, out_h), flags=flags, borderMode=cv2.BORDER_CONSTANT, borderValue=fill
    )
    if out.ndim == 2 and arr.ndim == 3:
        out = out[..., None]
    return out, {"mode": mode, "out_size": (out_w, out_h), "matrix": M.tolist(), "fill_value": list(fill)}


# --------------------------------------------------------------------------
# Text label (display copy only)
# --------------------------------------------------------------------------
COLORS = {
    "magenta": (255, 0, 255),
    "green": (0, 255, 0),
    "yellow": (255, 255, 0),
    "cyan": (0, 255, 255),
    "orange": (255, 128, 0),
    "red": (255, 0, 0),
    "blue": (0, 128, 255),
}


def parse_color(text: str):
    t = text.strip().lower()
    if t in COLORS:
        return COLORS[t]
    if t.startswith("#") and len(t) == 7:
        return tuple(int(t[i:i + 2], 16) for i in (1, 3, 5))
    parts = t.split(",")
    if len(parts) == 3:
        return tuple(max(0, min(255, int(x))) for x in parts)
    raise ValueError(f"bad color '{text}': use a name ({', '.join(COLORS)}), #RRGGBB, or R,G,B")


def to_display_rgb8(img: np.ndarray) -> np.ndarray:
    """Fixed linear conversion to 8-bit RGB (display only): uint16 -> /257, float -> min-max."""
    if img.dtype == np.uint8:
        a = img
    elif img.dtype == np.uint16:
        a = np.rint(img.astype(np.float32) / 257.0).astype(np.uint8)
    else:
        lo, hi = float(img.min()), float(img.max())
        a = np.rint((img - lo) / max(hi - lo, 1e-12) * 255).astype(np.uint8)
    if a.ndim == 2:
        a = np.stack([a, a, a], axis=-1)
    return np.ascontiguousarray(a)


DEFAULT_FONT = Path(__file__).resolve().parent / "fonts" / "lmsans10-bold.otf"  # Latin Modern Sans Bold


def load_font(path, size):
    from PIL import ImageFont
    return ImageFont.truetype(str(path), size)


def missing_glyphs(font, text):
    """Characters the font cannot draw (they would appear as blank/notdef boxes)."""
    sentinel = bytes(font.getmask("\ue000"))  # private-use char: absent from the font
    return sorted({c for c in text if not c.isspace() and bytes(font.getmask(c)) == sentinel})


def draw_label(img: np.ndarray, text: str, color, pad: int, font_path, size=None):
    """Return (labeled 8-bit RGB copy, info) with text at the gel's bottom-right."""
    from PIL import Image, ImageDraw

    rgb = to_display_rgb8(img)
    h, w = rgb.shape[:2]
    gel_w = max(w - 2 * pad, 1)
    size = int(size) if size else max(14, w // 28)
    while True:
        font = load_font(font_path, size)
        left, top, right, bottom = font.getbbox(text)
        if (right - left) <= 0.92 * gel_w or size <= 8:
            break
        size = max(8, int(size * 0.92))  # shrink long labels to fit inside the gel
    margin = max(6, int(0.01 * gel_w))
    pil = Image.fromarray(rgb)
    x0 = w - pad - margin - (right - left)   # right-align to the gel's bottom-right corner
    y0 = h - pad - margin - (bottom - top)
    ImageDraw.Draw(pil).text(
        (x0 - left, y0 - top), text, font=font,
        fill=tuple(int(c) for c in color),
    )
    info = {"font": " ".join(font.getname()), "font_file": Path(font_path).name,
            "font_size_px": size, "missing_glyphs": missing_glyphs(font, text)}
    return np.ascontiguousarray(np.asarray(pil)), info


# --------------------------------------------------------------------------
# Lane labels from comb templates (display copy only)
# --------------------------------------------------------------------------
COMBS_FILE = Path(__file__).resolve().parent / "combs.toml"


def load_combs(path=COMBS_FILE):
    with open(path, "rb") as f:
        return tomllib.load(f)


def comb_geometry(name, combs):
    """Return (well centres, pitch, well width), all as fractions of the gel width."""
    if name not in combs.get("combs", {}):
        raise ValueError(f"unknown comb '{name}'; available: {', '.join(combs.get('combs', {}))} "
                         f"(defined in {COMBS_FILE.name})")
    comb = combs["combs"][name]
    system = combs["systems"][comb["system"]]
    n = int(comb["wells"])
    ww = comb["well_width_mm"] / system["gel_width_mm"]
    first = system["first_well_left"] + ww / 2
    last = system["last_well_right"] - ww / 2
    if n == 1:
        return [(first + last) / 2], last - first, ww
    pitch = (last - first) / (n - 1)
    return [first + i * pitch for i in range(n)], pitch, ww


def parse_lanes(text: str, n_wells: int):
    """'8' -> ['1'..'8']; 'L,A1,,NTC' -> names (empty = unlabeled well), filled from the left."""
    t = text.strip()
    names = [str(i) for i in range(1, int(t) + 1)] if t.isdigit() else [s.strip() for s in t.split(",")]
    if len(names) > n_wells:
        raise ValueError(f"--lanes has {len(names)} entries but the comb has only {n_wells} wells")
    return names


def draw_lanes(rgb, names, centers_px, pitch_px, gel_top_px, font_path, size=None, color=None, angle="auto"):
    """Return (copy with lane labels in a new strip above the gel, info).

    angle: "auto" (horizontal if every label fits its lane, otherwise vertical), "0", "45" or "90".
    Rotated labels read upward, starting just above their well.
    """
    from PIL import Image, ImageDraw

    h, w = rgb.shape[:2]
    edge = np.concatenate([rgb[:2].reshape(-1, 3), rgb[:, :2].reshape(-1, 3), rgb[:, -2:].reshape(-1, 3)])
    bgc = tuple(int(v) for v in np.median(edge, axis=0))
    if color is None:
        color = (0, 0, 0) if sum(bgc) / 3 > 127 else (255, 255, 255)
    size = int(size) if size else max(10, int(pitch_px * 0.5))
    font = load_font(font_path, size)
    widths = [font.getbbox(n)[2] - font.getbbox(n)[0] if n else 0 for n in names]
    x0, y0, x1, y1 = font.getbbox("Ag0")
    text_h = max(y1 - y0, 1)
    if angle == "auto":
        theta = 90 if max(widths, default=0) > 0.9 * pitch_px else 0
    else:
        theta = int(angle)
    ct, st = math.cos(math.radians(theta)), math.sin(math.radians(theta))
    gap = max(4, size // 3)
    block = int(math.ceil(max(widths, default=0) * st + text_h * ct)) if theta else text_h
    top = block + 2 * gap
    right = 0
    if theta:
        right = max([int(math.ceil(x + wd * ct + text_h)) - w for x, wd in zip(centers_px, widths)] + [0])

    canvas = Image.new("RGB", (w + right, h + top), bgc)
    canvas.paste(Image.fromarray(rgb), (0, top))
    draw = ImageDraw.Draw(canvas)
    y_bottom = top + gel_top_px - gap
    fill = tuple(int(c) for c in color)
    for name, x, wd in zip(names, centers_px, widths):
        if not name:
            continue
        if not theta:
            draw.text((x, y_bottom), name, font=font, fill=fill, anchor="md")
            continue
        tw, th = wd + 2, text_h * 2
        tile = Image.new("RGBA", (tw, th), (0, 0, 0, 0))
        ImageDraw.Draw(tile).text((0, th / 2), name, font=font, fill=fill + (255,), anchor="lm")
        rot = tile.rotate(theta, resample=Image.BICUBIC, expand=True)
        px = rot.width / 2 - tw / 2 * ct   # where the start of the text ended up
        py = rot.height / 2 + tw / 2 * st
        canvas.paste(rot, (int(round(x - px)), int(round(y_bottom - py))), rot)
    info = {
        "font_size_px": size,
        "angle_deg": theta,
        "color": list(fill),
        "added_top_px": top,
        "added_right_px": right,
        "missing_glyphs": missing_glyphs(font, "".join(names)),
    }
    return np.ascontiguousarray(np.asarray(canvas)), info


# --------------------------------------------------------------------------
# QC image
# --------------------------------------------------------------------------
def _display(arr):
    a = arr
    if a.dtype != np.uint8:
        lo, hi = np.percentile(a, [0.5, 99.5])
        a = np.clip((a.astype(np.float32) - lo) / max(hi - lo, 1e-6) * 255, 0, 255).astype(np.uint8)
    if a.ndim == 2:
        a = cv2.cvtColor(a, cv2.COLOR_GRAY2BGR)
    else:
        a = cv2.cvtColor(a, cv2.COLOR_RGB2BGR)
    return a


def write_debug(path, arr, det, out, residual, lanes=None):
    """lanes: optional (centres_px, well_width_px, gel_top_px, gel_height_px) in output coordinates."""
    target_h = 700
    left = _display(arr).copy()
    s1 = target_h / left.shape[0]
    pts = (det["corners"] * s1).astype(np.int32).reshape(-1, 1, 2)
    left = cv2.resize(left, None, fx=s1, fy=s1, interpolation=cv2.INTER_AREA)
    cv2.polylines(left, [pts], True, (0, 200, 0), 2, cv2.LINE_AA)
    cv2.putText(left, f"detected tilt {det['phi']:+.2f} deg", (10, 28),
                cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 120, 255), 2, cv2.LINE_AA)
    right = _display(out)
    s2 = target_h / right.shape[0]
    right = cv2.resize(right, None, fx=s2, fy=s2, interpolation=cv2.INTER_AREA)
    if lanes:
        centres, ww, gtop, gh = lanes
        y_a, y_b = int(gtop * s2), int((gtop + 0.15 * gh) * s2)
        for x in centres:  # predicted well outlines from the comb template
            xa, xb = int((x - ww / 2) * s2), int((x + ww / 2) * s2)
            cv2.rectangle(right, (xa, y_a), (xb, y_b), (255, 0, 200), 1, cv2.LINE_AA)
    cv2.putText(right, f"output; residual tilt {residual:+.3f} deg", (10, 28),
                cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 120, 255), 2, cv2.LINE_AA)
    cv2.imwrite(str(path), np.hstack([left, np.full((target_h, 6, 3), 128, np.uint8), right]))


# --------------------------------------------------------------------------
# Per-file driver
# --------------------------------------------------------------------------
def process_file(src: Path, outdir: Path, args) -> str:
    stem = src.stem
    out_path = outdir / f"{stem}{args.suffix}.tif"
    log_path = outdir / f"{stem}{args.suffix}.json"
    if out_path.resolve() == src.resolve():
        raise ValueError("refusing to overwrite the input file")
    label_path = outdir / f"{stem}{args.suffix}_labeled.tif"
    preview_path = outdir / f"{stem}{args.suffix}_preview.png"
    want_labeled = bool(args.label or args.lanes)
    # --preview always reprocesses, so edited labels in a notebook block take effect.
    if (out_path.exists() and not args.force and not args.preview
            and (not want_labeled or label_path.exists())):
        return f"SKIP  {src.name} (output exists; use --force to redo)"
    combs = comb_info = None
    if args.lanes or args.comb:
        if not (args.lanes and args.comb):
            raise ValueError("--lanes and --comb must be used together")
        combs = load_combs()
        comb_info = comb_geometry(args.comb, combs)
        lane_names = parse_lanes(args.lanes, len(comb_info[0]))

    arr, res, notes = read_tiff(src)
    det = detect_gel(arr)
    out, info = rotate_crop(arr, det, args.pad, INTERP[args.interp], args.angle_tol)

    # Self-check: re-detect on the output; the tilt should now be ~0.
    try:
        residual = detect_gel(out)["phi"]
    except Exception:
        residual = float("nan")
    warnings = list(det["warnings"])
    if not math.isnan(residual) and abs(residual) > 0.25 and info["mode"] != "integer_crop_no_rotation":
        warnings.append(f"residual tilt after correction is {residual:+.2f} deg")

    outdir.mkdir(parents=True, exist_ok=True)
    provenance = {
        "tool": f"gel_autoalign {__version__}",
        "source_file": src.name,
        "rotation_deg": round(det["phi"], 4),
        "mode": info["mode"],
        "interpolation": args.interp,
    }
    write_tiff(out_path, out, res, json.dumps(provenance, separators=(",", ":")))

    # Gel position inside the output image (the slab is centred by rotate_crop).
    gel_w, gel_h = det["size"]
    out_w, out_h = info["out_size"]
    gel_left = (out_w - 1) / 2 - gel_w / 2
    gel_top = (out_h - 1) / 2 - gel_h / 2

    label_rec = lanes_rec = None
    lanes_dbg = None
    display = None
    if want_labeled:
        font_path = Path(args.font) if args.font else DEFAULT_FONT
        if not font_path.is_file():
            raise FileNotFoundError(f"label font not found: {font_path} (use --font, or keep the fonts/ folder next to the script)")
        display = to_display_rgb8(out)
        prov_l = dict(provenance, note="8-bit RGB display copy with burned-in text")
        text = None
        if args.label:
            text = args.label.replace("{name}", stem)
            display, finfo = draw_label(display, text, parse_color(args.label_color), args.pad, font_path, args.label_size)
            if finfo["missing_glyphs"]:
                warnings.append("label font has no glyph for: " + " ".join(finfo["missing_glyphs"]))
            prov_l["label"] = text
        if comb_info:
            fracs, pitch, ww = comb_info
            centres = [gel_left + f * gel_w for f in fracs]
            lanes_dbg = (centres, ww * gel_w, gel_top, gel_h)
            color = parse_color(args.lane_color) if args.lane_color else None
            display, linfo = draw_lanes(display, lane_names, centres[:len(lane_names)], pitch * gel_w,
                                        gel_top, font_path, args.lane_size, color, args.lane_angle)
            if linfo["missing_glyphs"]:
                warnings.append("lane font has no glyph for: " + " ".join(linfo["missing_glyphs"]))
            prov_l["lanes"] = lane_names
            comb = combs["combs"][args.comb]
            lanes_rec = {
                "comb": args.comb,
                "comb_spec": comb,
                "system_spec": combs["systems"][comb["system"]],
                "names": lane_names,
                "well_centres_px_in_aligned": [round(c, 2) for c in centres],
                "well_width_px": round(ww * gel_w, 2),
                **linfo,
                "note": "positions come from the comb template scaled to the detected gel width, not from the image",
            }
        write_tiff(label_path, display, res, json.dumps(prov_l, separators=(",", ":")))
        if args.label:
            label_rec = {
                "text": text,
                "color": args.label_color,
                **finfo,
            }
        file_rec = {"file": str(label_path.resolve()), "sha256": sha256(label_path),
                    "note": "display copy: 8-bit RGB, text burned in; the plain aligned TIFF is unchanged"}
        label_rec = {**label_rec, **file_rec} if label_rec else None
        if lanes_rec:
            lanes_rec.update(file_rec)

    preview_rec = None
    if args.preview:
        from PIL import Image
        img = Image.fromarray(display if display is not None else to_display_rgb8(out))
        pw = args.preview_width
        if img.width > pw:  # sized for a notebook window; the TIFFs keep full resolution
            img = img.resize((pw, round(img.height * pw / img.width)), Image.LANCZOS)
        img.save(preview_path)
        preview_rec = {"file": str(preview_path.resolve()), "note": "display-only PNG for notebooks"}

    record = {
        "tool": f"gel_autoalign {__version__}",
        "script_git": git_info(),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {
            "path": str(src.resolve()),
            "sha256": sha256(src),
            "shape": list(arr.shape),
            "dtype": str(arr.dtype),
        },
        "output": {"path": str(out_path.resolve()), "sha256": sha256(out_path), "shape": list(out.shape)},
        "parameters": {"pad_px": args.pad, "interpolation": args.interp, "angle_tol_deg": args.angle_tol},
        "detection": {
            "rotation_applied_deg_ccw": det["phi"],
            "min_area_rect_angle_deg": det["phi_minarearect"],
            "edge_fit_deviations_deg": det["side_deviations_deg"],
            "gel_center_px": det["center"],
            "gel_size_px": det["size"],
            "gel_corners_px": det["corners"].tolist(),
            "fill_ratio": det["fill_ratio"],
            "core_fill_ratio": det["core_fill_ratio"],
            "gel_area_fraction": det["area_fraction"],
            "background_gray": det["background_gray"],
        },
        "transform": info,
        "residual_tilt_after_deg": residual,
        "label": label_rec,
        "lanes": lanes_rec,
        "preview": preview_rec,
        "notes": notes,
        "warnings": warnings,
        "versions": {
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "opencv": cv2.__version__,
            "tifffile": tifffile.__version__,
        },
        "status": "warn" if warnings else "ok",
    }
    log_path.write_text(json.dumps(record, indent=2))
    if args.debug:
        write_debug(outdir / f"{stem}{args.suffix}_debug.png", arr, det, out, residual, lanes_dbg)

    tag = "WARN" if warnings else "OK  "
    line = (
        f"{tag}  {src.name}: rotated {det['phi']:+.2f} deg, "
        f"{arr.shape[1]}x{arr.shape[0]} -> {out.shape[1]}x{out.shape[0]}"
    )
    for wmsg in warnings:
        line += f"\n        ! {wmsg}"
    return line


def preview_link(path: Path) -> str:
    """Path for a notebook link: relative to the working directory (the org file's folder
    when run from org-babel) if the preview is inside it, absolute otherwise."""
    p = path.resolve()
    try:
        return str(p.relative_to(Path.cwd().resolve()))
    except ValueError:
        return str(p)


def collect(inputs, suffix):
    files = []
    for item in inputs:
        p = Path(item)
        if p.is_dir():
            files += sorted(
                f for f in p.iterdir()
                if f.is_file() and f.suffix.lower() in (".tif", ".tiff") and not f.stem.endswith(suffix)
            )
        elif p.is_file():
            files.append(p)
        else:
            print(f"? not found: {item}", file=sys.stderr)
    return files


def main(argv=None):
    ap = argparse.ArgumentParser(description="Automatically straighten and crop gel TIFFs.")
    ap.add_argument("inputs", nargs="+", help="TIFF file(s) and/or folder(s)")
    ap.add_argument("-o", "--outdir", help="output folder (default: 'aligned' next to each input)")
    ap.add_argument("--pad", type=int, default=20, help="margin around the gel in pixels (default 20)")
    ap.add_argument("--interp", choices=INTERP, default="cubic", help="interpolation for rotation (default cubic)")
    ap.add_argument("--angle-tol", type=float, default=0.05,
                    help="below this tilt (deg) do an exact integer crop with no resampling (default 0.05)")
    ap.add_argument("--suffix", default="_aligned", help="output filename suffix (default _aligned)")
    ap.add_argument("--label", help="text to print at the gel's bottom-right on an extra labeled copy; "
                    "'{name}' is replaced by the file name")
    ap.add_argument("--label-color", default="magenta",
                    help="label color: name, #RRGGBB or R,G,B (default magenta)")
    ap.add_argument("--label-size", type=int, default=None,
                    help="label font size in pixels (default: automatic, shrinks to fit)")
    ap.add_argument("--font", default=None,
                    help="font file for the label (default: fonts/lmsans10-bold.otf, Latin Modern Sans Bold)")
    ap.add_argument("--comb", help="comb template for lane labels (see combs.toml), e.g. mini-sub-15")
    ap.add_argument("--lanes", help="lane labels, filled from the left: a count ('8' -> 1..8) or names "
                    "separated by commas ('L,A1,,NTC'; an empty entry skips a well). Needs --comb")
    ap.add_argument("--lane-size", type=int, default=None,
                    help="lane label font size in pixels (default: automatic from the well spacing)")
    ap.add_argument("--lane-color", default=None,
                    help="lane label color: name, #RRGGBB or R,G,B (default: black or white, whichever contrasts)")
    ap.add_argument("--lane-angle", choices=["auto", "0", "45", "90"], default="auto",
                    help="lane label rotation: auto (horizontal if all labels fit, else vertical), "
                    "0, 45 or 90 degrees (default auto)")
    ap.add_argument("--preview", action="store_true",
                    help="also write a PNG preview and print only its path (for org-mode/notebooks); "
                    "always reprocesses")
    ap.add_argument("--preview-width", type=int, default=800,
                    help="maximum width of the --preview PNG in pixels (default 800)")
    ap.add_argument("--debug", action="store_true", help="also write a side-by-side QC PNG")
    ap.add_argument("--force", action="store_true", help="reprocess even if output exists")
    ap.add_argument("--watch", action="store_true", help="keep watching the input folder(s) for new files")
    args = ap.parse_args(argv)

    previews = []

    def run(files):
        n_warn = 0
        for f in files:
            outdir = Path(args.outdir) if args.outdir else f.parent / "aligned"
            try:
                msg = process_file(f, outdir, args)
            except Exception as exc:  # keep going on bad files
                msg = f"FAIL  {f.name}: {exc}"
                n_warn += 1
            n_warn += msg.startswith("WARN")
            # With --preview, stdout carries only the preview path (an org-babel file link), and
            # stderr stays empty unless there is a problem: org-babel pops up any stderr output.
            if not args.preview:
                print(msg, flush=True)
            elif not msg.startswith(("OK", "SKIP")):
                print(msg, flush=True, file=sys.stderr)
            preview = outdir / f"{f.stem}{args.suffix}_preview.png"
            if args.preview and not msg.startswith("FAIL") and preview.exists():
                previews.append(preview_link(preview))
        return n_warn

    def emit_previews():
        # No trailing newline: org-babel would put it inside the [[file:...]] link.
        if previews:
            sys.stdout.write("\n".join(previews))
            sys.stdout.flush()

    if not args.watch:
        files = collect(args.inputs, args.suffix)
        if not files:
            print("no TIFF files found", file=sys.stderr)
            return 1
        n_warn = run(files)
        emit_previews()
        return 0 if n_warn == 0 else 2

    print("watching for new TIFFs (Ctrl+C to stop)...", flush=True)
    sizes, done = {}, set()
    try:
        while True:
            for f in collect(args.inputs, args.suffix):
                if f in done:
                    continue
                size = f.stat().st_size
                if sizes.get(f) == size and size > 0:  # unchanged since last poll -> fully written
                    run([f])
                    done.add(f)
                sizes[f] = size
            time.sleep(2)
    except KeyboardInterrupt:
        print("stopped")
    return 0


if __name__ == "__main__":
    sys.exit(main())
