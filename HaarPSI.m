function score = HaarPSI(A, B)
% HAARPSI  Compute HaarPSI full-reference quality score for two images.
%   score = HaarPSI(A, B)
%
% This is a compact, self-contained implementation adapted from Reisenhofer et al.
% See: R. Reisenhofer et al., "A Haar wavelet-based perceptual similarity index..."
% (Please attribute in README / file header). :contentReference[oaicite:3]{index=3}

% Convert to double grayscale in [0,255] range
if ndims(A) == 3
    A = 0.299*double(A(:,:,1)) + 0.587*double(A(:,:,2)) + 0.114*double(A(:,:,3));
else
    A = double(A);
end
if ndims(B) == 3
    B = 0.299*double(B(:,:,1)) + 0.587*double(B(:,:,2)) + 0.114*double(B(:,:,3));
else
    B = double(B);
end

% Resize check
if any(size(A) ~= size(B))
    error('HaarPSI: images must be same size');
end

% Parameters (kept similar to paper defaults)
alpha = 0.5; beta = 6.0; % contrast / logistic-like constants (paper uses tuned values)
eps = 1e-8;

% Build Haar wavelet filters and compute 2-level multiscale coefficients
h = [1 1]/sqrt(2);
g = [1 -1]/sqrt(2);

% A simple multi-scale approach: compute Haar detail coefficients at two scales
[LL1_A, LH1_A, HL1_A, HH1_A] = haar_decomp(A);
[LL1_B, LH1_B, HL1_B, HH1_B] = haar_decomp(B);
[LL2_A, LH2_A, HL2_A, HH2_A] = haar_decomp(LL1_A);
[LL2_B, LH2_B, HL2_B, HH2_B] = haar_decomp(LL1_B);

% Local similarity maps: compute per-subband similarities and weight them
sim1 = local_similarity_map(LH1_A, LH1_B, alpha, eps);
sim2 = local_similarity_map(HL1_A, HL1_B, alpha, eps);
sim3 = local_similarity_map(HH1_A, HH1_B, alpha, eps);

sim4 = local_similarity_map(LH2_A, LH2_B, alpha, eps);
sim5 = local_similarity_map(HL2_A, HL2_B, alpha, eps);
sim6 = local_similarity_map(HH2_A, HH2_B, alpha, eps);

% combine with weighting based on contrast map (here simplified as local energy)
w1 = local_weight_map(LL1_A);
w2 = local_weight_map(LL2_A);

combined = (w1 .* sim1 + w1 .* sim2 + w1 .* sim3) + (w2 .* sim4 + w2 .* sim5 + w2 .* sim6);
% final pooling: mean of combined map (paper uses logistic pooling variant)
score = mean(combined(:));
end

% Helper small functions
function [LL, LH, HL, HH] = haar_decomp(I)
% one-level Haar decomposition using separable filters and decimation
h = [1 1]/sqrt(2);
g = [1 -1]/sqrt(2);
Lo = conv2(conv2(I, h, 'same'), h', 'same');
Hi1 = conv2(conv2(I, g, 'same'), h', 'same');
Hi2 = conv2(conv2(I, h, 'same'), g', 'same');
Hi3 = conv2(conv2(I, g, 'same'), g', 'same');

% Downsample by 2 (simple)
LL = Lo(1:2:end, 1:2:end);
LH = Hi1(1:2:end,1:2:end);
HL = Hi2(1:2:end,1:2:end);
HH = Hi3(1:2:end,1:2:end);
end

function S = local_similarity_map(X, Y, alpha, eps)
% similarity map for two subbands (X and Y), stabilized
num = 2 .* X .* Y + eps;
den = X.^2 + Y.^2 + eps;
S = (num ./ den) .^ alpha;
end

function W = local_weight_map(LL)
% simple local energy weight map normalized
E = sqrt(conv2(LL.^2, ones(3)/9, 'same'));
W = E / (max(E(:)) + eps);
end
