import numpy as np


def _center_crop(A, marg=0.1):
    """Remove a 10% border to avoid boundary artifacts (mirrors MATLAB centerCrop)."""
    h, w = A.shape
    mh = max(1, round(h * marg))
    mw = max(1, round(w * marg))
    # MATLAB A(mh:h-mh+1, ...) is 1-indexed → Python A[mh-1:h-mh+1, ...]
    return A[mh - 1:h - mh + 1, mw - 1:w - mw + 1]


def dlmaim(R_coeffs, O_coeffs, A_coeffs):
    """
    Detail Loss Measure with Additive Impairment Modulation.

    Parameters
    ----------
    R_coeffs : list — abs of restored coefficients,  pywt.wavedec2 structure
    O_coeffs : list — abs of reference coefficients, pywt.wavedec2 structure
    A_coeffs : list — abs of additive-impairment coefficients, same structure

    Each list: [approx_arr, (H_n, V_n, D_n), ..., (H_1, V_1, D_1)]

    Returns
    -------
    q : float — DLMAIM quality score
    """
    eps = np.finfo(np.float64).eps

    # ------------------------------------------------------------------
    # Part 1: detail restoration term  (L3-norm ratio across all subbands)
    # ------------------------------------------------------------------
    nume1 = 0.0
    deno1 = 0.0
    for i in range(1, len(R_coeffs)):        # skip approximation band (index 0)
        for theta in range(3):               # H=0, V=1, D=2
            Rsub = _center_crop(R_coeffs[i][theta])
            Osub = _center_crop(O_coeffs[i][theta])
            nume1 += float(np.sum(Rsub.ravel() ** 3) ** (1.0 / 3.0))
            deno1 += float(np.sum(Osub.ravel() ** 3) ** (1.0 / 3.0))

    q1 = nume1 / (deno1 + eps)

    # ------------------------------------------------------------------
    # Part 2: additive impairment penalty
    # ------------------------------------------------------------------
    # Total coefficient count (approximation + all detail subbands)
    total_size = A_coeffs[0].size + sum(
        A_coeffs[i][theta].size
        for i in range(1, len(A_coeffs))
        for theta in range(3)
    )

    nume2 = 0.0
    for i in range(1, len(A_coeffs)):
        for theta in range(3):
            Asub = _center_crop(A_coeffs[i][theta])
            SD = float(np.std(Asub))
            nume2 += float(np.sum((Asub.ravel() + SD) ** 3) ** (1.0 / 3.0))

    a1 = 10.0
    q2 = nume2 / (total_size + eps)

    return q1 - a1 * q2
