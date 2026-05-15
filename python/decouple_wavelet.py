import numpy as np


def _cal_angle(H, V):
    return np.arctan(V / (H + 1e-30)) + np.pi * (H < 0)


def decouple_wavelet(coeffs_ref, coeffs_test):
    """
    Wavelet-domain decoupling of detail losses and additive impairments.

    Parameters
    ----------
    coeffs_ref  : list — pywt.wavedec2 output for the reference image
    coeffs_test : list — pywt.wavedec2 output for the test/distorted image

    Returns
    -------
    R : list — restored (detail-loss) coefficients, same structure as input
    A : list — additive-impairment coefficients  (= coeffs_test - R)
    """
    # Copy approximation band from test image
    R = [coeffs_test[0].copy()]

    for i in range(1, len(coeffs_ref)):
        OH, OV, OD = coeffs_ref[i]
        TH, TV, TD = coeffs_test[i]

        # Gain factor clamped to [0, 1] — models detail loss
        k_H = np.clip(TH / (OH + 1e-30), 0.0, 1.0)
        k_V = np.clip(TV / (OV + 1e-30), 0.0, 1.0)
        k_D = np.clip(TD / (OD + 1e-30), 0.0, 1.0)

        RH = k_H * OH
        RV = k_V * OV
        RD = k_D * OD

        # Contrast-enhancement correction for H and V subbands:
        # where the gradient orientation is preserved, use test coefficients directly.
        OA = _cal_angle(OH, OV)
        TA = _cal_angle(TH, TV)
        pos = np.abs(OA - TA) * (180.0 / np.pi) < 1.0

        RH[pos] = TH[pos]
        RV[pos] = TV[pos]

        R.append((RH, RV, RD))

    # Additive impairment is whatever the restored signal didn't account for
    A = [coeffs_test[0] - R[0]]
    for i in range(1, len(coeffs_test)):
        TH, TV, TD = coeffs_test[i]
        RH, RV, RD = R[i]
        A.append((TH - RH, TV - RV, TD - RD))

    return R, A
