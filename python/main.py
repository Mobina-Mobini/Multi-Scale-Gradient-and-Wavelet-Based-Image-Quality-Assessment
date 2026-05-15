import numpy as np
import pywt

from gr import gr
from decouple_wavelet import decouple_wavelet
from dlmaim import dlmaim
# from haarpsi import haarpsi  # available but disabled by default (see haarpsi.py)


def _rgb2gray(img):
    img = np.asarray(img, dtype=np.float64)
    return 0.299 * img[:, :, 0] + 0.587 * img[:, :, 1] + 0.114 * img[:, :, 2]


def _abs_coeffs(coeffs):
    """Return a new coefficient list with absolute values applied throughout."""
    result = [np.abs(coeffs[0])]
    for i in range(1, len(coeffs)):
        result.append(tuple(np.abs(sub) for sub in coeffs[i]))
    return result


def main(ref_img, test_img):
    """
    Compute the MSGW quality score between a reference and a distorted image.

    Parameters
    ----------
    ref_img  : numpy array — reference image (H×W grayscale or H×W×3 RGB)
    test_img : numpy array — distorted image, same shape as ref_img

    Returns
    -------
    score : float — combined MSGW + DLMAIM quality score
                    (higher means the test image is closer to the reference)

    Example
    -------
    >>> from PIL import Image
    >>> import numpy as np
    >>> ref  = np.array(Image.open('ref.png'))
    >>> dist = np.array(Image.open('distorted.png'))
    >>> score = main(ref, dist)
    >>> print(score)
    """
    # Grayscale conversion
    ref_y  = _rgb2gray(ref_img)  if ref_img.ndim  == 3 else np.asarray(ref_img,  dtype=np.float64)
    test_y = _rgb2gray(test_img) if test_img.ndim == 3 else np.asarray(test_img, dtype=np.float64)

    if ref_y.shape != test_y.shape:
        raise ValueError("Reference and test images must have the same size.")

    levels = 4
    wname  = 'db2'
    mode   = 'periodization'   # equivalent to MATLAB dwtmode('per')

    # Wavelet decomposition
    coeffs_r = pywt.wavedec2(ref_y,  wname, level=levels, mode=mode)
    coeffs_t = pywt.wavedec2(test_y, wname, level=levels, mode=mode)

    # Multi-scale gradient similarity
    # gr() is computed at the original scale and at each of the 4 approximation levels
    sim_product = gr(ref_y, test_y)
    for lev in range(1, levels + 1):
        Ar = pywt.wavedec2(ref_y,  wname, level=lev, mode=mode)[0]
        At = pywt.wavedec2(test_y, wname, level=lev, mode=mode)[0]
        sim_product *= gr(Ar, At)

    # Wavelet decoupling → DLMAIM
    R_coeffs, A_coeffs = decouple_wavelet(coeffs_r, coeffs_t)
    dlmaim_score = dlmaim(
        _abs_coeffs(R_coeffs),
        _abs_coeffs(coeffs_r),
        _abs_coeffs(A_coeffs),
    )

    # HaarPSI (optional — enable below if needed)
    # try:
    #     haar_score = haarpsi(ref_y, test_y)
    # except Exception as e:
    #     haar_score = 1.0
    #     print(f"HaarPSI skipped: {e}")

    return sim_product * dlmaim_score
    # return sim_product * dlmaim_score * haar_score


# ---------------------------------------------------------------------------
# Quick smoke-test when run directly:  python main.py [ref] [distorted]
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import sys
    from PIL import Image

    if len(sys.argv) == 3:
        ref  = np.array(Image.open(sys.argv[1]))
        test = np.array(Image.open(sys.argv[2]))
        print(f"Score: {main(ref, test):.6f}")
    else:
        rng   = np.random.default_rng(0)
        img   = rng.integers(0, 256, (256, 256, 3), dtype=np.uint8)
        noisy = np.clip(img.astype(np.int32) + rng.integers(-20, 20, img.shape), 0, 255).astype(np.uint8)
        print(f"Same image  score: {main(img, img):.6f}")
        print(f"Noisy image score: {main(img, noisy):.6f}")
