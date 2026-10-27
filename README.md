# Fixed-emitter displacement estimation

MATLAB implementation accompanying *Computational Measurement of Nanoscale Displacement from Fixed-Emitter Image Sequences*.

The algorithm detects and localizes emitters, matches them to a reference template, rejects outliers using MAD, and estimates displacement. Template updates and coordinate bridging provide a global trajectory.

## Requirements

Tested with MATLAB R2020b, Image Processing Toolbox, and Statistics and Machine Learning Toolbox. A C++ compiler is optional.

## Quick start

For a self-contained demonstration, open this directory in MATLAB and run:

```matlab
[results, truth, trajectoryFig] = example_run;
```

This generates 25 small TIFFs in a new temporary directory and displays the estimated trajectory. The fixed-seed example uses Gaussian emitters and additive noise; it demonstrates the workflow rather than the paper's camera-noise model. The image directory is printed in the Command Window.

For your own images:

1. Open `main_drift_beads.m` in MATLAB.
2. Set `params.dataDir` to a TIFF sequence directory and set `params.pixelSizeNm` to its calibration.
3. Adjust the detection, matching, and template parameters in the script, then click **Run**.

Keep `private/` beside the entry files; do not add it separately to the MATLAB path. Settings for your data are defined in the main script; the example has its own fixed settings. Use the figure-specific settings when reproducing paper results.

Input images must be same-size, single-channel TIFFs. Use zero-padded filenames in frame order, such as `frame_0001.tif`; filenames are sorted alphabetically. At least five frames with sufficient stable emitters are needed for initialization. Manual detection thresholds use normalized image intensities (`uint16` values divided by 65535).

The [image datasets](../DATA_README.md) are stored separately. For Fig. 2, select an individual `frames/` directory.

## Output

The script displays a 2D global trajectory with template-update markers and leaves `params`, `results`, and `trajectoryFig` in the workspace.

| Result field | Description |
| --- | --- |
| `traj_global` | Global displacement in nm, relative to frame 1 |
| `step_pair_nm` | Displacement between consecutive frames in nm |
| `rel_to_template_px` | Displacement relative to the active template in pixels |
| `frame_n_used` | Number of matches retained after MAD rejection in each frame |
| `frame_match_ok` | Measurement success flags; failed frames hold the previous position |
| `updateFrames`, `update_events` | Successful template switches and details of update attempts |

Coordinates follow image columns (x) and rows (y). Displacement is current minus reference; an image correction uses the opposite sign.

```matlab
save('displacement_results.mat', 'results', 'params');
exportgraphics(trajectoryFig, 'trajectory.png', 'Resolution', 300);
```

## Optional MEX acceleration

The code uses the MATLAB localizer when the package-local MEX binary is unavailable. To compile it, run the following from this directory with a supported C++ compiler:

```matlab
privateDir = fullfile(pwd, 'private');
mex('-R2018a', '-outdir', privateDir, ...
    fullfile(privateDir, 'localize_weighted_centroid_mex.cpp'));
```
