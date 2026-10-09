"""
Unified quality metrics for the cryo-ET reconstruction/denoising benchmark.

Shared by all three tracks (denoise / reconstruct / pick). Pure NumPy/SciPy/
scikit-image so it runs on Windows with no GPU.

Metrics
-------
- snr_db / snr_regions : global and region-based signal-to-noise ratio
- cnr                  : contrast-to-noise ratio between two regions
- psnr, ssim           : full-reference (need ground truth)
- fsc, fsc_resolution  : Fourier shell correlation + resolution at a threshold
- pick_prf             : precision / recall / F1 for 3D particle picking

All functions are documented with the exact definition used, so the numbers in
the paper's Table 2 are reproducible.
"""
from __future__ import annotations
import numpy as np


# ----------------------------------------------------------------------------- SNR / CNR
def snr_db(signal: np.ndarray, noise_std: float) -> float:
    """20*log10(RMS(signal)/noise_std)  [dB]."""
    rms = float(np.sqrt(np.mean(np.square(signal))))
    return 20.0 * np.log10(rms / noise_std) if noise_std > 0 else np.inf


def snr_regions(vol: np.ndarray, signal_mask: np.ndarray,
                background_mask: np.ndarray) -> float:
    """Region SNR = (mean_signal - mean_bg) / std_bg.

    signal_mask / background_mask are boolean arrays of vol's shape.
    """
    sig = vol[signal_mask]
    bg = vol[background_mask]
    denom = bg.std()
    return float((sig.mean() - bg.mean()) / denom) if denom > 0 else np.inf


def cnr(vol: np.ndarray, mask_a: np.ndarray, mask_b: np.ndarray) -> float:
    """Contrast-to-noise ratio = |mean_a - mean_b| / sqrt(var_a + var_b)."""
    a, b = vol[mask_a], vol[mask_b]
    denom = np.sqrt(a.var() + b.var())
    return float(abs(a.mean() - b.mean()) / denom) if denom > 0 else np.inf


# ----------------------------------------------------------------------------- edge sharpness
def edge_sharpness(vol: np.ndarray) -> float:
    """Mean Sobel gradient magnitude over the volume.

    Sensitive to blur: smoothing lowers it even when region-SNR rises, which is
    exactly the failure mode of naive SNR scoring (a Gaussian-blurred volume
    scores high SNR but low sharpness). Computed on the z-scored volume so the
    value is comparable across differently scaled outputs.
    """
    from scipy.ndimage import sobel
    v = (vol - vol.mean()) / (vol.std() + 1e-8)
    g2 = np.zeros_like(v, dtype=np.float64)
    for ax in range(v.ndim):
        d = sobel(v, axis=ax)
        g2 += d.astype(np.float64) ** 2
    return float(np.sqrt(g2).mean())


# ----------------------------------------------------------------------------- full reference
def psnr(pred: np.ndarray, gt: np.ndarray, data_range: float | None = None) -> float:
    from skimage.metrics import peak_signal_noise_ratio
    if data_range is None:
        data_range = float(gt.max() - gt.min())
    return float(peak_signal_noise_ratio(gt, pred, data_range=data_range))


def ssim(pred: np.ndarray, gt: np.ndarray, data_range: float | None = None) -> float:
    from skimage.metrics import structural_similarity
    if data_range is None:
        data_range = float(gt.max() - gt.min())
    return float(structural_similarity(gt, pred, data_range=data_range))


# ----------------------------------------------------------------------------- FSC
def fsc(vol_a: np.ndarray, vol_b: np.ndarray, n_shells: int | None = None,
        voxel_size=1.0, power_floor=1e-12):
    """FSC on physical-frequency shells, including non-cubic arrays.

    Each axis uses fftfreq(n, d); voxel_size is a scalar or (z,y,x).
    Frequencies are cycles per spacing unit, limited to the smallest axial
    Nyquist frequency. DC is excluded. rFFT coefficients have Hermitian
    multiplicities, so sums equal a full FFT. Empty/negligible-power shells
    are NaN, not artificial zero correlation. No masking or taper is added.
    """
    from scipy.fft import rfftn, fftfreq, rfftfreq
    a, b = np.asarray(vol_a), np.asarray(vol_b)
    if a.shape != b.shape or a.ndim != 3 or min(a.shape) < 4:
        raise ValueError('FSC requires equal three-dimensional shapes, axes >=4')
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('FSC inputs must be finite')
    spacing = np.broadcast_to(np.asarray(voxel_size, float), (3,))
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError('voxel_size must be finite and positive')
    n_shells = min(a.shape)//2 if n_shells is None else int(n_shells)
    if n_shells < 2:
        raise ValueError('At least two frequency shells are required')
    nyquist = float(np.min(.5/spacing))
    width = nyquist/n_shells
    freq = (np.arange(n_shells)+.5)*width
    # SciPy preserves single precision, avoiding huge complex128 allocations.
    fa, fb = rfftn(a, workers=2), rfftn(b, workers=2)
    fz, fy = fftfreq(a.shape[0],spacing[0]), fftfreq(a.shape[1],spacing[1])
    fx = rfftfreq(a.shape[2],spacing[2])
    yz = fy[:,None]**2+fx[None,:]**2
    mult = np.full(fx.size,2.0); mult[0]=1
    if a.shape[2]%2==0: mult[-1]=1
    weight=np.broadcast_to(mult,yz.shape)
    num=np.zeros(n_shells); pa=num.copy(); pb=num.copy(); count=num.copy()
    for z in range(a.shape[0]):
        radius=np.sqrt(fz[z]**2+yz)
        valid=(radius>0)&(radius<=nyquist)
        idx=np.minimum((radius[valid]/width).astype(int),n_shells-1)
        az,bz=fa[z][valid],fb[z][valid]; w=weight[valid]
        num+=np.bincount(idx,weights=(az*np.conj(bz)).real*w,minlength=n_shells)
        pa+=np.bincount(idx,weights=np.abs(az)**2*w,minlength=n_shells)
        pb+=np.bincount(idx,weights=np.abs(bz)**2*w,minlength=n_shells)
        count+=np.bincount(idx,weights=w,minlength=n_shells)
    curve=np.full(n_shells,np.nan)
    ok=(count>0)&(pa>pa.max()*power_floor)&(pb>pb.max()*power_floor)
    curve[ok]=np.clip(num[ok]/np.sqrt(pa[ok]*pb[ok]),-1,1)
    return freq,curve


def fsc_crossing(freq, curve, threshold=0.143, pixel_size_a=None):
    """First bracketed downward crossing; never substitute the last bin.

    no_crossing: correlation remains above threshold on the valid contiguous
    band, frequency >= frequency_bound, resolution <= resolution_bound_A.
    below_first_shell: unresolved low-frequency crossing; no point estimate.
    invalid_gap/no_valid_shells: cannot infer a crossing through absent power.
    pixel_size_a converts cycles/pixel to Angstrom; omit for physical grids.
    """
    freq,curve=np.asarray(freq,float),np.asarray(curve,float)
    if freq.shape!=curve.shape or np.any(np.diff(freq)<=0):
        raise ValueError('Frequencies must be strictly increasing and match curve')
    ids=np.flatnonzero(np.isfinite(curve)&np.isfinite(freq)&(freq>0))
    result=dict(status='no_valid_shells',frequency=np.nan,frequency_bound=np.nan,
                resolution_A=np.nan,resolution_bound_A=np.nan,threshold=float(threshold))
    if not len(ids): return result
    first=int(ids[0]); last=first
    if curve[first]<threshold:
        result.update(status='below_first_shell',frequency_bound=float(freq[first]))
    else:
        for i in range(first+1,len(freq)):
            if not np.isfinite(curve[i]): break
            if curve[i]<threshold:
                f=freq[i-1]+(curve[i-1]-threshold)/(curve[i-1]-curve[i])*(freq[i]-freq[i-1])
                result.update(status='crossing',frequency=float(f))
                break
            last=i
        else: last=len(freq)-1
        if result['status']!='crossing':
            result.update(status='invalid_gap' if np.any(ids>last) else 'no_crossing',
                          frequency_bound=float(freq[last]))
    if pixel_size_a is not None:
        if np.isfinite(result['frequency']): result['resolution_A']=float(pixel_size_a/result['frequency'])
        if np.isfinite(result['frequency_bound']): result['resolution_bound_A']=float(pixel_size_a/result['frequency_bound'])
    return result


def fsc_resolution(vol_a, vol_b, threshold: float = 0.143, pixel_size_a: float | None = None):
    """Bracketed crossing in cycles/pixel or Angstrom; NaN when unbracketed.

    Use fsc_crossing for censoring status and bounds. Legacy endpoint-valued
    results produced before 2026-09-16 must be recomputed, not rescaled.
    """
    info=fsc_crossing(*fsc(vol_a,vol_b),threshold=threshold,pixel_size_a=pixel_size_a)
    return info['resolution_A'] if pixel_size_a is not None else info['frequency']


# ----------------------------------------------------------------------------- particle picking
def pick_prf(pred_xyz: np.ndarray, gt_xyz: np.ndarray, tol: float):
    """Greedy one-to-one matching within Euclidean distance `tol` (pixels/nm).

    Returns dict with precision, recall, f1, tp, fp, fn, mean_loc_error.
    """
    pred_xyz = np.asarray(pred_xyz, float).reshape(-1, 3)
    gt_xyz = np.asarray(gt_xyz, float).reshape(-1, 3)
    if len(pred_xyz) == 0 or len(gt_xyz) == 0:
        tp = 0
        errs = []
    else:
        from scipy.spatial import cKDTree
        tree = cKDTree(gt_xyz)
        used = set()
        tp = 0
        errs = []
        # sort predictions arbitrarily; greedy nearest assignment
        for p in pred_xyz:
            d, j = tree.query(p, k=min(5, len(gt_xyz)))
            d = np.atleast_1d(d); j = np.atleast_1d(j)
            for dj, jj in zip(d, j):
                if dj <= tol and jj not in used:
                    used.add(int(jj)); tp += 1; errs.append(float(dj)); break
    fp = len(pred_xyz) - tp
    fn = len(gt_xyz) - tp
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0
    return dict(precision=precision, recall=recall, f1=f1, tp=tp, fp=fp, fn=fn,
                mean_loc_error=float(np.mean(errs)) if errs else np.nan)


# ----------------------------------------------------------------------------- self-test
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    # FSC sanity: identical volumes -> ~1 everywhere; independent noise -> ~0
    v = rng.standard_normal((48, 48, 48))
    f, c = fsc(v, v)
    print(f"[self-test] FSC(v,v) mean = {c.mean():.3f} (expect ~1.0)")
    f, c = fsc(v, rng.standard_normal((48, 48, 48)))
    print(f"[self-test] FSC(v,noise) mean = {c.mean():.3f} (expect ~0.0)")
    # PSNR/SSIM sanity
    gt = rng.standard_normal((64, 64))
    noisy = gt + 0.5 * rng.standard_normal((64, 64))
    print(f"[self-test] PSNR(noisy) = {psnr(noisy, gt):.2f} dB, SSIM = {ssim(noisy, gt):.3f}")
    # picking sanity: perfect prediction -> P=R=F1=1
    gtp = rng.uniform(0, 100, (20, 3))
    print(f"[self-test] pick perfect: {pick_prf(gtp, gtp, tol=1.0)}")
    # picking with half missing -> recall 0.5
    print(f"[self-test] pick half:    {pick_prf(gtp[:10], gtp, tol=1.0)}")
    # edge sharpness sanity on a STRUCTURED volume (the intended use case):
    # blurring a signal-bearing volume must lower it. (On pure noise the metric
    # is dominated by the noise itself and is not meaningful in isolation --
    # that is exactly why it must be read jointly with SNR, see paper Sec. 4.)
    from scipy.ndimage import gaussian_filter
    zz, yy, xx = np.mgrid[0:48, 0:48, 0:48]
    ball = (((zz - 24) ** 2 + (yy - 24) ** 2 + (xx - 24) ** 2) <= 12 ** 2).astype(float)
    ball += 0.05 * rng.standard_normal(ball.shape)
    es_raw = edge_sharpness(ball)
    es_blur = edge_sharpness(gaussian_filter(ball, 2.0))
    print(f"[self-test] edge_sharpness sphere={es_raw:.3f} blurred={es_blur:.3f} (expect sphere > blurred)")
    assert es_raw > es_blur, "edge_sharpness must decrease when a structured volume is blurred"
    print("[self-test] OK")
