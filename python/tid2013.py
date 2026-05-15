"""
TID2013 dataset evaluation using the MSGW quality metric.

Expected dataset layout:
    <data_path>/
        reference_images/   I01.BMP ... I25.BMP
        distorted_images/   I01_01_1.BMP ... I25_24_5.BMP
        mos.txt             (optional) plain-text MOS vector, one value per line

Usage:
    python tid2013.py /path/to/tid2013
"""

import os
import sys
import numpy as np
from PIL import Image
from scipy.stats import spearmanr

from main import main as msgw_score


def load_image(path):
    return np.array(Image.open(path))


def evaluate(data_path):
    ref_folder  = os.path.join(data_path, 'reference_images')
    dist_folder = os.path.join(data_path, 'distorted_images')
    mos_file    = os.path.join(data_path, 'mos.txt')

    if not os.path.isdir(ref_folder) or not os.path.isdir(dist_folder):
        raise FileNotFoundError(
            f"Expected sub-folders not found:\n  {ref_folder}\n  {dist_folder}"
        )

    scores = np.zeros(3000)
    mld    = np.zeros((24, 25 * 5))
    i_point = 0

    for i_ref in range(1, 26):
        ref_path = os.path.join(ref_folder, f'I{i_ref:02d}.BMP')
        s = load_image(ref_path).astype(np.float64)
        lum_r = 0.299 * s[:, :, 0] + 0.587 * s[:, :, 1] + 0.114 * s[:, :, 2]

        for i_dis in range(1, 25):
            for i_level in range(1, 6):
                dist_path = os.path.join(
                    dist_folder, f'I{i_ref:02d}_{i_dis:02d}_{i_level}.BMP'
                )
                if not os.path.isfile(dist_path):
                    raise FileNotFoundError(f"Distorted image missing: {dist_path}")

                t = load_image(dist_path).astype(np.float64)
                lum_t = 0.299 * t[:, :, 0] + 0.587 * t[:, :, 1] + 0.114 * t[:, :, 2]

                lum_diff = float(np.max(lum_t - lum_r))
                mld[i_dis - 1, (i_ref - 1) * 5 + (i_level - 1)] = lum_diff

                # Chroma-like Q component for moderate luminance shifts
                if 0.9 < lum_diff < 1.4:
                    q_r = 0.211 * s[:, :, 0] - 0.523 * s[:, :, 1] + 0.312 * s[:, :, 2]
                    q_t = 0.211 * t[:, :, 0] - 0.523 * t[:, :, 1] + 0.312 * t[:, :, 2]
                    scores[i_point] = msgw_score(q_r, q_t)
                else:
                    scores[i_point] = msgw_score(lum_r, lum_t)

                i_point += 1
                print(f"  [{i_point:4d}/3000] ref={i_ref:02d} dis={i_dis:02d} lv={i_level}  score={scores[i_point-1]:.5f}")

    results = {'scores': scores, 'mld': mld}

    if os.path.isfile(mos_file):
        mos = np.loadtxt(mos_file)
        rho, _ = spearmanr(scores, mos)
        print(f"\nOverall Spearman correlation: {rho:.4f}")
        results['mos'] = mos
        results['spearman'] = rho

        scores_mat = scores.reshape(25 * 5, 24)
        mos_mat    = mos.reshape(25 * 5, 24)
        cor_adm = np.array([
            spearmanr(mos_mat[:, j], scores_mat[:, j])[0]
            for j in range(24)
        ])
        results['spearman_per_distortion'] = cor_adm
        print("Per-distortion Spearman:")
        print(np.round(cor_adm, 4))
    else:
        print(f"MOS file not found at {mos_file}. Raw scores saved only.")

    np.save('tid2013_results.npy', results)
    print("Results saved to tid2013_results.npy")
    return results


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else input("Enter TID2013 dataset path: ").strip()
    evaluate(path)
