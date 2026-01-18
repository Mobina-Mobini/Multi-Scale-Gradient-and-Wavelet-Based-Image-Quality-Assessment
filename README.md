# Multi-Scale-Gradient-and-Wavelet-Based-Image-Quality-Assessment
MATLAB implementation of an objective image quality assessment (IQA) method based on Multi-scale Gradient similarity and Wavelet.


The Multi-scale Gradient Wavelet (MSGW) method is a full-reference image quality assessment metric proposed by Mobini and Faraji. MSGW applies a 4-level Haar wavelet decomposition to the reference and distorted images and computes gradient-magnitude similarity at each scale. This captures both low-frequency (approximation) and high-frequency image information. A color-saturation distortion detection step is included, and the final MSGW score combines the multi-scale gradient–wavelet measure with the Haar Perceptual Similarity Index (HaarPSI) to better reflect human perception.

Code Structure

main.m: Main function ADMScore = main(refImg, distImg) that computes the MSGW score. It converts color images to grayscale, performs 4-level discrete Haar wavelet decomposition (wavedec2), and computes gradient similarity at each scale via the helper function gr. It also calls sub-functions (decouple, DLMAIM, etc.) to measure Depletion (detail-loss) and pure dissonance (additive impairments), and combines these with the HaarPSI score.

gr.m: Computes a gradient-similarity score between two images. It filters each image with a Scharr gradient kernel and computes a similarity index (analogous to GMSD) by pooling the lowest 20% of similarity values.

HaarPSI.m: Implements the Haar Wavelet-Based Perceptual Similarity Index as described by Reisenhofer et al.. This is called in main.m to obtain HaarScore.

Additional routines (in main.m):

decouple(O,T,Z): Splits the wavelet coefficients into a “detail loss” image R and an “additive impairment” image A.

get_subband_se, get_subband_2D, cal_angle, center, DLMAIM: Support functions that index and reshape subbands, compute angular contrast, crop central regions, and compute the detail-loss measure (DLMAIM) according to Eq. (19) of [24].

Data script: TID2013YIQ.m (included for preprocessing TID2013 dataset images into YIQ color space, if needed).

Usage

Place all files in the MATLAB path. Ensure the Wavelet Toolbox is available.

Call the main function with the reference and distorted images. For example:

ref = imread('ref.png');
dist = imread('distorted.png');
score = main(ref, dist);


Citation

If you use this code or the MSGW method, please cite the original paper: M. Mobini and M. R. Faraji, “Multi-scale Gradient Wavelet-based Image Quality Assessment,” Visual Computer, vol. 40, pp. 8713–8728, 2024.


