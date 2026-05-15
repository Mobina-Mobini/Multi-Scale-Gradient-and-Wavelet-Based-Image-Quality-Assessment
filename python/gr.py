import numpy as np
from scipy.signal import convolve2d


def gr(ref, test):
    """
    Gradient similarity pooling between two grayscale images.
    Returns a scalar in roughly [0, 1] (higher = more similar).
    """
    A = np.asarray(ref, dtype=np.float64)
    B = np.asarray(test, dtype=np.float64)

    # Scharr-like kernels (normalized)
    dx = np.array([[3, 0, -3],
                   [10, 0, -10],
                   [3, 0, -3]], dtype=np.float64) / 16.0
    dy = dx.T

    IxA = convolve2d(A, dx, mode='same')
    IyA = convolve2d(A, dy, mode='same')
    gA = np.sqrt(IxA ** 2 + IyA ** 2)

    IxB = convolve2d(B, dx, mode='same')
    IyB = convolve2d(B, dy, mode='same')
    gB = np.sqrt(IxB ** 2 + IyB ** 2)

    # Stabilized pointwise similarity
    epsilon = 1e-6
    T = 1000.0
    sim_map = (2.0 * gA * gB + T) / (gA ** 2 + gB ** 2 + T + epsilon)

    # Pool the worst 20% of similarity values
    r = 0.20
    n = max(1, round(r * sim_map.size))
    vals = np.sort(sim_map.ravel())
    return float(np.mean(vals[:n]))
