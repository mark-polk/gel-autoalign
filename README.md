# gel-autoalign

`gel_autoalign.py` automatically straightens and crops gel images (TIFF). It finds the gel slab with classical OpenCV image processing (blur → Otsu threshold → largest object → edge line fits), measures its tilt, and rotates and crops in a single resampling step. If the tilt is negligible, it does an exact integer crop with no interpolation.

It's deterministic, with no machine learning and no random numbers: the same input, parameters and library versions give a byte-identical output.

## Setup

```bash
mamba env create -f environment.yml
mamba activate gel-autoalign
python gel_autoalign.py --help
```

`environment.lock.yml` records the exact package versions the tool was last set up with. To reproduce that environment exactly, run `mamba env create -f environment.lock.yml`.

## Usage

```bash
python gel_autoalign.py gel.tif                     # one file
python gel_autoalign.py /path/to/folder             # every .tif/.tiff in the folder
python gel_autoalign.py /path/to/folder --watch     # keep processing new files as they appear
python gel_autoalign.py gel.tif --pad 40            # margin around the gel in pixels (default 20)
python gel_autoalign.py gel.tif --debug             # also write a side-by-side QC PNG
python gel_autoalign.py gel.tif --label "MP-02 PCR" # also write a labeled copy
python gel_autoalign.py folder --label "{name}"     # label each gel with its file name
```

Label options:

- `--label TEXT` prints the text at the gel's bottom-right on an extra copy. `{name}` is replaced by each input file's name (without extension).
- `--label-color` sets the label color as a name (`magenta`, the default, or `green`, `yellow`, `cyan`, `orange`, `red`, `blue`), as `#RRGGBB`, or as `R,G,B`.
- `--label-size` sets the font size in pixels. By default it's chosen automatically and shrunk to fit inside the gel.
- `--font` sets the font file to use. The default is `fonts/lmsans10-bold.otf` (Latin Modern Sans Bold, bundled).

Other options: `-o/--outdir`, `--interp {cubic,linear,lanczos}`, `--angle-tol`, `--suffix`, `--force`. Run `--help` for details.

## Outputs

Outputs are written to an `aligned/` folder next to each input, unless you set `-o`:

| File | Contents |
|---|---|
| `<name>_aligned.tif` | Straightened, cropped gel. Lossless (zlib), with the same bit depth, channels and DPI as the input. |
| `<name>_aligned.json` | Record of the run: input/output SHA-256, parameters, detected angle and geometry, transform matrix, warnings, library versions, and the script's git state (`script_git`). |
| `<name>_aligned_labeled.tif` | Only written with `--label`. A labeled display copy. |
| `<name>_aligned_debug.png` | Only written with `--debug`. Side-by-side QC image showing the detection and the result. |

## Caveats

- **Rotation resamples pixels.** For densitometry or any quantification, use the raw imager files, not the aligned output.
- **The script removes tilt but doesn't decide orientation.** It corrects tilt modulo 90° and doesn't know which side is "up", so place gels on the imager in a consistent orientation.
- **The labeled copy is an 8-bit RGB display file** with the text burned in. The plain `_aligned.tif` is never modified by labeling and remains the clean record.
- **Originals are never modified.** Outputs always go to a separate folder.
