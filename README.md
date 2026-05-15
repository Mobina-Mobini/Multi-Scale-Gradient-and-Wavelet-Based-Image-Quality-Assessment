# Multi-Scale Gradient Wavelet-Based Image Quality Assessment (MSGW)

A full-reference Image Quality Assessment (IQA) metric based on Multi-Scale Gradient similarity and Wavelet decomposition, proposed by **Mobini and Faraji**.

> **Citation:** M. Mobini and M. R. Faraji, "Multi-scale Gradient Wavelet-based Image Quality Assessment," *The Visual Computer*, vol. 40, pp. 8713-8728, 2024.

---

## Method Overview

MSGW applies a 4-level wavelet decomposition to reference and distorted images and computes gradient-magnitude similarity at each scale. This captures both low-frequency (approximation) and high-frequency detail information. The final score combines:

- **Multi-scale gradient similarity** -- Scharr-filtered gradient maps pooled over the worst 20% of local similarity values, repeated at each wavelet approximation level.
- **DLMAIM** -- Detail Loss Measure with Additive Impairment Modulation, which separates distortions into detail-loss and additive-noise components in the wavelet domain.
- **HaarPSI** *(optional)* -- Haar Wavelet-Based Perceptual Similarity Index (Reisenhofer et al.), available but disabled by default.

---

## Repository Structure

```
├── matlab/          # Original MATLAB implementation
│   ├── main.m
│   ├── gr.m
│   ├── decouple_wavelet.m
│   ├── dlmaim.m
│   ├── HaarPSI.m
│   └── TID2013YIQ.m
│
├── python/          # Python implementation (NumPy / PyWavelets)
│   ├── main.py
│   ├── gr.py
│   ├── decouple_wavelet.py
│   ├── dlmaim.py
│   ├── haarpsi.py
│   ├── tid2013.py
│   └── requirements.txt
│
└── README.md
```

---

## MATLAB Usage

Requires the **Wavelet Toolbox**. Place all files from `matlab/` on the MATLAB path.

```matlab
ref  = imread('ref.png');
dist = imread('distorted.png');
score = main(ref, dist);
```

---

## Python Usage

### Install dependencies

```bash
pip install -r python/requirements.txt
```

Dependencies: `numpy`, `scipy`, `PyWavelets`, `Pillow`

### Run on two images

```python
from PIL import Image
import numpy as np
import sys
sys.path.insert(0, 'python')
from main import main

ref   = np.array(Image.open('ref.png'))
dist  = np.array(Image.open('distorted.png'))
score = main(ref, dist)
print(score)
```

Or directly from the command line:

```bash
cd python
python main.py ref.png distorted.png
```

### Quick smoke-test (no images needed)

```bash
cd python
python main.py
```

### Evaluate on TID2013

```bash
cd python
python tid2013.py /path/to/tid2013
```

Expected dataset layout:
```
tid2013/
    reference_images/   I01.BMP ... I25.BMP
    distorted_images/   I01_01_1.BMP ... I25_24_5.BMP
    mos.txt             (optional)
```

---

## Notes on the Python Translation

- `pywt.wavedec2` with `mode='periodization'` replaces MATLAB's `wavedec2` + `dwtmode('per')`.
- The Python version works directly with PyWavelets' structured coefficient list, eliminating the need for MATLAB's flat coefficient vector and bookkeeping matrix.
- The HaarPSI size-inconsistency bug present in the original MATLAB code is fixed in `haarpsi.py` by upsampling level-2 maps before combining. It remains disabled by default in `main.py` to match the original behaviour.
