PENG EXTENSION, 2026-10-03
This is a bounded supplement to a review article, not a new algorithm or a
comprehensive matched-information method benchmark. No original result was replaced.

LAYOUT
The archive preserves experiments/revision_20261003_extension and the small
calibration dependencies in experiments/data, experiments/code and
experiments/revision_20261003_peng_followup. Extract into a NEW project folder.
The separate Micron_LaTeX_source.zip supplies the updated manuscript project.

RECHECK SAVED HALF RESULTS (no downloads or retraining)
Use a Python environment with numpy, scipy, pandas, mrcfile and scikit-image.
Run from the extension directory:
  python -X utf8 run_halves.py evaluate
This reads the bundled 96-cube input/prediction arrays and reproduces 24 rows
and 48 FSC curves. The full-FFT equivalence check must pass. Gaussian-tapered
spectra are calculated in float64 to avoid single-precision tail error.
Running validate_delivery.py additionally requires the two full manuscript
versions and the local cryoCARE training-data records; it is a local delivery
audit, not the portable half-array recheck command.

NEW DATA AND FROZEN CZII EXPERIMENTS
acquire.py fetches all five unused DS10440 runs and hybrid-curated point labels.
acquire.py halves fetches the exact Figshare tutorial archive. Downloaded data
must match retained hashes before using it as the same input version; portal
URLs/annotation inventories may change. Discovery metadata were captured in
portal_inventory.json. Download records use real retrieval timestamps, not
file modification dates. About 1.46 GB of new CZII volumes and a 1.95 GB tutorial
ZIP were downloaded. Raw full-volume data and environments are NOT bundled.
Use the recorded original environment/package lists and a compatible GPU for
Topaz; the local main environment uses Python 3.11 and PyTorch 2.6 CUDA12.4.
Run run_samples.py preprocess, then run_samples.py picking. They write only
the extension folder. Frozen thresholds come from the included previous
TS_5_4 calibration CSV; no new-run outcomes choose them. Existing caches are
reused, so use a NEW extension directory for an independent rerun.
The small-radius candidate-pruning change preserves selected peaks at the
frozen operating point; verify_peak_pruning.py records an exact coordinate
and ordering comparison with a cached original 0.30-floor extraction.

OFFICIAL CRYOCARE TRAINING
Actual package versions: cryocare0.2.2, csbdeep0.7.4, TensorFlow2.15.1,
NumPy1.26.4, Python3.11, CPU. Complete dependency list is retained.
The executed run_cryocare.py uses ASCII work directory
C:/Users/wangz/cryoet_extension_runtime_20261003 to avoid historical Windows
TensorFlow path issues. On another computer, change only its WORK path to a
NEW ASCII folder before running, record that path-only difference, and retain
all model/input settings. The existing completion files suppress retraining;
for an independent run use a fresh output directory, not the archived one.
The official package source is not modified. The wrapper chooses train/test
crops and calls its DataModule, CryoCARE.train and predict_direct.
Training crop: z[0,208), y[0,768), x[0,544); train y[0,651), validation
y[651,768), guard y[768,896). Fixed geometric test centers are in the protocol.
20 epochs x20 steps, batch2, patch64, depth3, first16channels, seed20261003.
The validation-selected checkpoint is epoch20. Training 310.551 s; inference
times are in cryocare_complete.json. This is NOT a convergence claim.
Weights/configuration and the sampled training/validation coordinates are
bundled in cryocare_model and cryocare_training_data; paths within the saved
NPZ files record the actual workstation and need relocation for dataset reload.
The pretrained DeepDeWedge tutorial checkpoint was NOT used.

EXECUTION RECORD / KNOWN ATTEMPTS
- An initial HEAD request to Figshare returned400; subsequent GET succeeded.
- Download connections timed out and retried; final lengths and hashes saved.
- First picking launch failed on the Windows default GBK text decoder before
  scoring. It was rerun with -X utf8, without changing scientific settings.
- The broad 0.30-floor small-radius peak extraction was interrupted to prune
  unneeded candidates BELOW the already frozen score thresholds. Existing
  candidate caches were preserved; no performance-selected parameter changed.
- Initial float32 Hann-tapered spectral tails failed the strict independent
  FFT check. Float64 resolved the numerical mismatch; tolerance was not relaxed.
- All final experiments and evaluations finished. No monitor or compute task
  was left running by this workflow.

INTERPRETATION LIMITS
Five CZII runs add acquisitions, not demonstrated biological preparations.
Their annotations were hybrid-curated using multiple processed inputs.
Tomo110 halves are frame-disjoint observations with common acquisition,
alignment and reconstruction; six spatial crops are ONE specimen. The MRC
voxel spacing is zero: no Angstrom resolution is inferred. Cross-half MSE is
prediction of the opposite noisy observation, not clean-reference MSE.
Possible Topaz pretraining overlap with Tomo110 is unknown. The fixed-budget
cryoCARE run does not isolate architecture, domain match or computational cost.
The Hann-window check was added after initial unwindowed curves were inspected;
both versions and the amendment record are retained. No gold-standard STA,
before-splitting leakage test or additional supervised cryo-ET method was run.

No public release or permanent DOI was created. Consult upstream data/code
licenses before redistribution; this local audit bundle is not a license grant.
