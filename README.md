# gel-autoalign

Automatically straightens and crops gel images (TIFF files).

Point it at a gel image, or a whole folder of them, and it will:

1. find the gel in the picture,
2. measure how crooked it is,
3. rotate it straight and crop away the empty background, and
4. save the result as a new file, leaving your original untouched.

It can also stamp a label (such as the file name) in the corner of the gel.

It uses standard image processing (OpenCV), not AI. The same image always gives exactly the same result.

> **A friendly heads-up:** this was vibe coded (written with an AI assistant) for quick internal use, such as tidying up gels for lab notebooks, lab meetings and day-to-day sharing. It works well for that, but it hasn't been rigorously validated. Always look over the results yourself.
>
> When preparing gel images to share outside the lab, for a paper, a poster or a collaborator, take a bit more care. Keep the raw imager files, say what processing you did, and check the image guidelines of the journal or venue you're submitting to.

---

## One-time setup

You only need to do this once per computer. It takes about 10 minutes.

### Step 1: Install Miniforge

Miniforge is a free installer that sets up Python and the tools this script needs, without touching anything else on your computer.

**Mac**

1. Open the **Terminal** app. Press ⌘+Space, type `Terminal`, then press Enter.
2. Copy and paste this line into Terminal, then press Enter:

   ```bash
   curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
   ```

3. Then paste this line and press Enter:

   ```bash
   bash Miniforge3-$(uname)-$(uname -m).sh
   ```

4. Press Enter to scroll through the license, type `yes` to accept, and press Enter to accept the default install location. When it asks whether to initialize ("conda init"), type `yes`.
5. **Close the Terminal window and open a new one.** Your prompt should now start with `(base)`.

**Windows**

1. Download the Windows installer (`Miniforge3-Windows-x86_64.exe`) from the [Miniforge download page](https://github.com/conda-forge/miniforge#download).
2. Run it and accept the defaults.
3. From the Start menu, open **Miniforge Prompt**. On Windows, type every command in this guide into the Miniforge Prompt, not the regular Command Prompt.

**Linux:** follow the Mac steps in your usual terminal.

> Already have Anaconda, Miniconda or Miniforge? Skip to Step 2.

### Step 2: Download this tool

1. On this GitHub page, click the green **Code** button, then **Download ZIP**.
2. Unzip it and move the `gel-autoalign` folder somewhere you'll remember, such as your Documents folder.

(If you use git, you can run `git clone` instead.)

### Step 3: Install the tool's software

1. In Terminal (or Miniforge Prompt), type `cd` followed by a space. **Don't press Enter yet.**
2. Drag the `gel-autoalign` folder from Finder (or File Explorer) into the Terminal window. Its location appears after `cd`. Now press Enter.
3. Paste this line and press Enter:

   ```bash
   conda env create -f environment.yml
   ```

   This downloads everything the script needs into a self-contained environment called `gel-autoalign`. It takes a few minutes. When the prompt comes back, setup is done.

---

## Using it

### 1. Open Terminal and turn on the environment

You need to do this each time you open a new Terminal window:

```bash
conda activate gel-autoalign
```

Your prompt now starts with `(gel-autoalign)` instead of `(base)`.

### 2. Run the script on your gels

The easiest way is to drag things into the Terminal window instead of typing their locations:

1. Type `python` followed by a space.
2. Drag the file **`gel_autoalign.py`** (inside the `gel-autoalign` folder) into the window, then type a space.
3. Drag in **your gel image**, or **a folder of gel images**.
4. Press Enter.

The line will look something like this:

```bash
python /Users/you/Documents/gel-autoalign/gel_autoalign.py /Users/you/Desktop/my_gels
```

For each image you'll see a line like:

```
OK    gel1.tif: rotated +1.37 deg, 2048x1536 -> 1702x1210
```

`OK` means it worked. `WARN` means it worked but something looked unusual, so check that image by eye. `FAIL` means it couldn't process that image. The reason is printed next to it.

### 3. Find your results

The results go into a new folder called **`aligned`**, created next to your original images. **Your original files are never changed.**

## Common options

Add these at the end of the command, after the image or folder:

| What you want | Add this |
|---|---|
| Label each gel with its own file name | `--label "{name}"` |
| Label with your own text | `--label "MP-02 PCR"` |
| Change the label color | `--label-color yellow` (also `magenta` (the default), `green`, `cyan`, `orange`, `red`, `blue`, or a hex code like `#FF8800`) |
| Make the label bigger or smaller | `--label-size 60` (in pixels; by default it's sized automatically) |
| Use a different font | `--font /path/to/font.ttf` |
| Leave more margin around the gel | `--pad 40` (in pixels; the default is 20) |
| Label the lanes (see [Lane labels](#lane-labels)) | `--comb mini-sub-15 --lanes "1kb,A1,A2,NTC"` |
| Save a PNG preview and print its path (for notebooks) | `--preview` |
| Also save a "before / after" check image | `--debug` |
| Keep running and process new images as they appear in a folder | `--watch` (press Ctrl+C to stop) |
| Redo images that were already processed | `--force` |

You can combine them. For example, this processes every gel in a folder, labels each one with its file name, and saves check images. (As before, you can drag in the script and the folder instead of typing them.)

```bash
python gel_autoalign.py my_gels --label "{name}" --debug
```

The quotes around `"{name}"` are needed.

To see every option, run `python gel_autoalign.py --help`.

## Lane labels

Tell the script which comb you used and what's in each lane. It works out where the wells are from the comb's layout and writes the names above the gel.

```bash
python gel_autoalign.py gel.tif --comb mini-sub-15 --lanes "1kb,WT,ΔabcA,,NTC"
```

- **`--lanes`** lists the names from left to right, separated by commas. Leave an entry empty (`,,`) to skip an unloaded well. Wells after the last name stay unlabeled. Typing just a number, such as `--lanes 8`, labels the lanes 1 to 8.
- Short names are written horizontally. If any name is too wide for its lane, all the labels are written at 45°.
- **`--lane-size`** sets the font size. **`--lane-color`** sets the color; by default it's black or white, whichever stands out against the background.

Combs included:

| `--comb` | Comb | Gel system |
|---|---|---|
| `mini-protean-10` | Bio-Rad 10-well, 1.0 mm (#1653359) | Mini-PROTEAN (PAGE) |
| `mini-protean-15` | Bio-Rad 15-well, 1.0 mm (#1653360) | Mini-PROTEAN (PAGE) |
| `mini-sub-8` | Bio-Rad 8-well, 1.5 mm (#1704463) | Mini-Sub Cell GT (agarose) |
| `mini-sub-15` | Bio-Rad 15-well, 1.5 mm (#1704465) | Mini-Sub Cell GT (agarose) |

**How it works.** The well positions are measured as a fraction of the gel's width, so they don't depend on the camera or zoom. The templates were measured from real gels made with the 10-well PAGE and 15-well agarose combs. The other two combs are calculated from those, assuming every comb in a system starts and ends its wells at the same place.

**Check your first gel with each comb.** Run it with `--debug`. The check image outlines where the script expects each well in magenta, so you can see whether they line up with the real wells. This matters most for `mini-protean-15` and `mini-sub-8`, which haven't been checked against a real gel yet. If a gel was trimmed or damaged at the sides, the positions will shift with it.

To add another comb, copy an entry in `combs.toml` and edit it. The comments in that file explain each value.

## Using it from Emacs org-mode

`--preview` saves a PNG of the result and prints only its path. That's the format org-babel needs to insert the image into your notebook:

```org
#+begin_src sh :results file link
conda run -n gel-autoalign python ~/tools/gel-autoalign/gel_autoalign.py \
  ~/notebook/gels/2026-11-13_PCR.tif --comb mini-sub-15 --lanes "1kb,A1,A2,NTC" --preview
#+end_src
```

Press `C-c C-c` in the block. To see the image itself rather than a link, turn on inline images with `M-x org-display-inline-images`, or add `#+STARTUP: inlineimages` to the file. Any messages, such as warnings, are printed separately so they don't break the link.

With `--preview`, the image is always reprocessed, so after editing the lane names you can just run the block again. The full-quality TIFFs and the `.json` record are saved to `aligned/` as usual. Point the command at a copy of the imager file, not the only original.

## What you get

For an input called `gel1.tif`, the `aligned` folder will contain:

| File | What it is |
|---|---|
| `gel1_aligned.tif` | The straightened, cropped gel. It keeps the full quality of the original (same bit depth and resolution, lossless compression). |
| `gel1_aligned.json` | A record of exactly what was done: the rotation angle, settings, file checksums and software versions. Keep it with your data. |
| `gel1_aligned_labeled.tif` | Only made if you used `--label` or `--lanes`. A labeled copy for slides and figures. |
| `gel1_aligned_preview.png` | Only made if you used `--preview`. A display-only PNG for notebooks. |
| `gel1_aligned_debug.png` | Only made if you used `--debug`. A side-by-side image showing what the script detected and the result. |

## Important limitations

- **Don't use the output for quantification.** Rotating an image slightly changes its pixel values. For densitometry or any band quantification, always use the raw files from the imager.
- **It straightens, but it doesn't flip.** It removes tilt, but it can't tell which side of the gel is "up". Place gels on the imager in the same orientation every time.
- **The labeled copy is for display only.** It's converted to an 8-bit color image with the text burned in. The plain `_aligned.tif` is never labeled and stays the clean record.
- **Originals are never modified.** Results always go into a separate folder.
- **The `.json` record contains the full folder paths of your files.** Keep that in mind before sharing it.

## Troubleshooting

| Problem | Fix |
|---|---|
| `conda: command not found` | Close Terminal and open a new window. If it still happens, re-run the Miniforge installer and answer `yes` to the "conda init" question. On Windows, make sure you're in **Miniforge Prompt**. |
| `Shell not initialized` after typing `mamba activate` | Use `conda activate gel-autoalign` instead. |
| `No module named cv2` (or `numpy`, `tifffile`, `PIL`) | You forgot to turn on the environment. Run `conda activate gel-autoalign`. |
| `can't open file 'gel_autoalign.py'` | Terminal can't find the script. Drag `gel_autoalign.py` into the window instead of typing its name. |
| `SKIP ... output exists` | That image was already processed. Add `--force` to redo it. |
| A gel came out cropped wrong | Run it again with `--debug` and look at the `_debug.png` image. The script works best when the gel clearly stands out from an even background and doesn't touch the edge of the photo. |

## Updating

Download the new ZIP (or run `git pull`), then turn on the environment and update it:

```bash
conda activate gel-autoalign
conda env update -f environment.yml
```

---

## For the technically inclined

- The script is a single file, `gel_autoalign.py`. How it detects the gel: Gaussian blur → Otsu threshold → largest object → robust line fits to the four edges, starting from OpenCV's minimum-area rectangle. Rotation and cropping happen in one resampling step (cubic by default; set with `--interp`). If the tilt is below `--angle-tol`, it does an exact integer crop with no interpolation.
- `environment.yml` is the file to install from, on any operating system. `environment.lock.yml` records the exact versions last tested on macOS (Apple Silicon). It's there for reproducibility and won't install on other platforms.
- Lane positions come from the templates in `combs.toml`. For each gel system it stores where the first well's left edge and the last well's right edge fall, as fractions of the gel width. Wells are spaced evenly between those two edges, and each well's width is taken from the manufacturer's spec. The computed positions and the template used are written to the `.json` record. The script only places labels; it doesn't detect lanes in the image.
- When the script is run from a git checkout, each `.json` record includes the script's git commit under `script_git`.
- `.gitignore` excludes all TIFFs and `aligned/` folders, so gel images can't be committed by accident.

## License and credits

The code is released under the MIT License (see `LICENSE`). You're free to use, modify and share it.

The label font is Latin Modern Sans Bold, by GUST e-foundry. It's included in `fonts/` under the GUST Font License (see `fonts/GUST-FONT-LICENSE.TXT`).
