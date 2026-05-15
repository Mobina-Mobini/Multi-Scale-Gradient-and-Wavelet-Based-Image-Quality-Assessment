import numpy as np
from scipy.signal import convolve2d


def _haar_decomp(I):
    """One-level Haar decomposition with downsampling by 2."""
    h = np.array([[1.0, 1.0]]) / np.sqrt(2)
    g = np.array([[1.0, -1.0]]) / np.sqrt(2)

    Lo  = convolve2d(convolve2d(I, h,  mode='same'), h.T, mode='same')
    Hi1 = convolve2d(convolve2d(I, g,  mode='same'), h.T, mode='same')
    Hi2 = convolve2d(convolve2d(I, h,  mode='same'), g.T, mode='same')
    Hi3 = convolve2d(convolve2d(I, g,  mode='same'), g.T, mode='same')

    return Lo[::2, ::2], Hi1[::2, ::2], Hi2[::2, ::2], Hi3[::2, ::2]


def _local_similarity_map(X, Y, alpha=0.5, eps=1e-8):
    num = 2.0 * X * Y + eps
    den = X ** 2 + Y ** 2 + eps
    return (num / den) ** alpha


def _local_weight_map(LL, eps=1e-8):
    E = np.sqrt(convolve2d(LL ** 2, np.ones((3, 3)) / 9.0, mode='same'))
    return E / (E.max() + eps)


def _upsample(x, target_shape):
    """Nearest-neighbour upsample to target_shape."""
    rh = max(1, target_shape[0] // x.shape[0])
    rw = max(1, target_shape[1] // x.shape[1])
    return np.repeat(np.repeat(x, rh, axis=0), rw, axis=1)[:target_shape[0], :target_shape[1]]


def haarpsi(A, B):
    """
    Haar Wavelet-Based Perceptual Similarity Index.
    Adapted from Reisenhofer et al.

    The original MATLAB implementation had a size inconsistency when combining
    level-1 and level-2 similarity maps. This Python version fixes that by
    upsampling level-2 maps to level-1 size before pooling.

    Parameters
    ----------
    A, B : 2-D or 3-D numpy arrays (grayscale or RGB, any numeric dtype)

    Returns
    -------
    score : float
    """
    def to_gray(img):
        img = np.asarray(img, dtype=np.float64)
        if img.ndim == 3:
            return 0.299 * img[:, :, 0] + 0.587 * img[:, :, 1] + 0.114 * img[:, :, 2]
        return img

    A = to_gray(A)
    B = to_gray(B)

    if A.shape != B.shape:
        raise ValueError("HaarPSI: images must be the same size.")

    alpha = 0.5
    eps = 1e-8

    LL1_A, LH1_A, HL1_A, HH1_A = _haar_decomp(A)
    LL1_B, LH1_B, HL1_B, HH1_B = _haar_decomp(B)
    LL2_A, LH2_A, HL2_A, HH2_A = _haar_decomp(LL1_A)
    LL2_B, LH2_B, HL2_B, HH2_B = _haar_decomp(LL1_B)

    sim1 = _local_similarity_map(LH1_A, LH1_B, alpha, eps)
    sim2 = _local_similarity_map(HL1_A, HL1_B, alpha, eps)
    sim3 = _local_similarity_map(HH1_A, HH1_B, alpha, eps)

    sim4 = _local_similarity_map(LH2_A, LH2_B, alpha, eps)
    sim5 = _local_similarity_map(HL2_A, HL2_B, alpha, eps)
    sim6 = _local_similarity_map(HH2_A, HH2_B, alpha, eps)

    w1 = _local_weight_map(LL1_A)
    w2 = _local_weight_map(LL2_A)

    target = sim1.shape
    combined = (
        w1 * sim1 + w1 * sim2 + w1 * sim3
        + _upsample(w2, target) * _upsample(sim4, target)
        + _upsample(w2, target) * _upsample(sim5, target)
        + _upsample(w2, target) * _upsample(sim6, target)
    )

    return float(np.mean(combined))
