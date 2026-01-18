function s = gr(ref, test)
% GR Gradient similarity pooling between two (grayscale) images.
%   s = gr(ref, test)
%
% Uses Scharr-like kernels to obtain gradient magnitude maps, computes
% pointwise similarity map and pools using the mean of the worst r fraction.
%
% Returns scalar s in [0,1] roughly (higher = more similar).

% Ensure double
A = double(ref);
B = double(test);

% Scharr-like filter kernels (normalized)
dx = [3 0 -3; 10 0 -10; 3 0 -3] / 16;
dy = dx';

% gradient maps
IxA = conv2(A, dx, 'same');
IyA = conv2(A, dy, 'same');
gA = sqrt(IxA.^2 + IyA.^2);

IxB = conv2(B, dx, 'same');
IyB = conv2(B, dy, 'same');
gB = sqrt(IxB.^2 + IyB.^2);

% pointwise stabilized similarity
epsilon = 1e-6;            % small stabilization
T = 1000;                  % large stabilizer (keeps values in dynamic range)
simMap = (2 .* gA .* gB + T) ./ (gA.^2 + gB.^2 + T + epsilon);

% pooling: mean of lowest r fraction (emphasize worst local gradient similarity)
r = 0.20;
n = max(1, round(r * numel(simMap)));
vals = sort(simMap(:), 'ascend');
s = mean(vals(1:n));
end
