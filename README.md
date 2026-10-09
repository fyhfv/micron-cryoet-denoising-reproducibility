# Reproducibility package — cryo-ET denoising / missing-wedge review (Micron)

Bounded **supplementary reproducibility materials** for the review manuscript:

> *Deep-Learning Denoising and Missing-Wedge Restoration in Cryo-Electron Tomography: A Review of Methods, Benchmarks and Downstream Consequences*

This archive supports the **Peng_extension** controls (2026-10-03). It is **not** a new restoration algorithm and **not** a comprehensive matched-information method benchmark. No original result in the review was replaced by these controls.

> **Zenodo DOI:** *not yet assigned* — after the first GitHub Release is linked to Zenodo, replace this line with the version DOI (e.g. `https://doi.org/10.5281/zenodo.XXXXXXX`) and update the manuscript *Data and code availability* section.

---

## What you can recheck without downloads or retraining

From `experiments/revision_20261003_extension/`:

```bash
python -X utf8 run_halves.py evaluate
```

Requirements: Python with `numpy`, `scipy`, `pandas`, `mrcfile`, `scikit-image`.

This reads the bundled half-set cubes, reproduces **24 metric rows** and **48 FSC curves**, and must pass the full-FFT equivalence check (Gaussian-tapered spectra in float64).

`validate_delivery.py` is a **local delivery audit** (needs full manuscript versions and local cryoCARE training-data records). It is not the portable half-array recheck.

---

## Layout

| Path | Contents |
|------|----------|
| `experiments/revision_20261003_extension/` | Extension scripts, half evaluation arrays, cryoCARE weights/logs, provenance, manifests |
| `experiments/revision_20261003_peng_followup/` | Calibration / follow-up dependencies used by the extension |
| `experiments/code/`, `experiments/data/` | Shared metric code and small calibration metadata |
| `summary_tables/` | Convenience copies of key CSV/JSON tables used in the manuscript package |
| `README_reproducibility_original.txt` | Original delivery README (verbatim) |

Extract / clone into a **new** working folder. Do not overwrite an existing experiment tree if you need an independent rerun.

---

## New CZII runs and frozen thresholds (optional, downloads required)

`acquire.py` fetches five unused DS10440 runs and hybrid-curated point labels; `acquire.py halves` fetches the Figshare DeepDeWedge tutorial archive. Downloaded bytes must match retained hashes before treating them as the same input version.

- Portal dataset: [CZII dataset 10440](https://cryoetdataportal.czscience.com/datasets/10440)
- Tutorial archive: [Figshare file 45582309](https://api.figshare.com/v2/file/download/45582309)

Raw full-volume data and conda/venv environments are **not** fully bundled. Recorded environments:

- Main / Topaz-oriented: Python 3.11, PyTorch 2.6 CUDA 12.4 (see `environment_main.json`)
- Official cryoCARE fit: cryocare 0.2.2, csbdeep 0.7.4, TensorFlow 2.15.1, NumPy 1.26.4, CPU (`environment_cryocare.json`)

Frozen picking thresholds come from the included previous TS_5_4 calibration CSV; **new-run outcomes must not choose them**.

---

## Official cryoCARE demonstration

The wrapper calls the unmodified official package (`DataModule`, `CryoCARE.train`, `predict_direct`) with fixed crops and seed `20261003`. See `cryocare_protocol.json` and `run_cryocare.py`. On another machine, set only the ASCII `WORK` path to a **new** folder; keep model/input settings unchanged. Existing completion files suppress retraining—use a fresh output directory for an independent run.

---

## Interpretation limits (please keep)

- Five CZII runs add acquisitions, not demonstrated independent biological preparations.
- Tomo110 halves are frame-disjoint observations with shared acquisition / alignment / reconstruction.
- No Ångström resolution is inferred from these half-data (MRC voxel spacing recorded as zero in the tutorial arrays).
- No gold-standard experimental STA, before-splitting training-leakage control, or permanent claim of method ranking is provided by this package alone.

---

## Upstream data and licenses

Public tomograms and tutorial reconstructions remain under their **original** portal / Figshare / software licenses. This repository redistributes **scripts, hashes, derived tables, selected arrays, and fitted demonstration weights** for audit and recheck. It is **not** a license grant for bulk re-hosting of all upstream raw volumes. Consult CZII, Figshare, cryoCARE, DeepDeWedge, and Topaz terms before further redistribution.

Code authored for this extension package is released under the MIT License (see `LICENSE`). Third-party code and data keep their original terms (see `NOTICE`).

---

## Citation

Until a Zenodo DOI exists, cite the manuscript and this repository URL. After Zenodo archival, prefer the **version DOI** of the release that matches the submitted/accepted text.

A `CITATION.cff` file is included for GitHub’s citation panel.

---

## Maintainer notes (authors)

1. Connect this repository in [Zenodo → Settings → GitHub](https://zenodo.org/account/settings/github/).
2. Enable the repo toggle **before** creating a GitHub Release.
3. Create Release tag `v1.0.0` (or similar).
4. Paste the Zenodo version DOI into this README and into the manuscript *Data and code availability* paragraph (replace “accompanying reproducibility files” / “identifier has not yet been assigned”).
5. For double-blind review, use a restricted Zenodo record or anonymous link if the journal requires it; do not invent a DOI.
