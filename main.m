function score = main(refImg, testImg)
% MAIN Compute combined MSGW + ADM(DLMAIM) + HaarPSI quality score

% -------------------------
% Input handling
% -------------------------
if ndims(refImg) == 3
    refY = rgb2gray_uint8safe(refImg);
else
    refY = double(refImg);
end

if ndims(testImg) == 3
    testY = rgb2gray_uint8safe(testImg);
else
    testY = double(testImg);
end

if any(size(refY) ~= size(testY))
    error('Reference and test images must have the same size.');
end

% -------------------------
% Parameters
% -------------------------
levels = 4;
wname  = 'db2';

% -------------------------
% Wavelet decomposition
% -------------------------
dwtmode('per','nodisp');
[coeffR, bookR] = wavedec2(refY, levels, wname);
[coeffT, bookT] = wavedec2(testY, levels, wname);

% -------------------------
% Multi-scale gradient similarity
% -------------------------
simProduct = gr(refY, testY);
for lev = 1:levels
    Ar = appcoef2(coeffR, bookR, wname, lev);
    At = appcoef2(coeffT, bookT, wname, lev);
    simProduct = simProduct * gr(Ar, At);
end

% -------------------------
% Decoupling + DLMAIM
% -------------------------
[Rcoeff, Acoeff] = decouple_wavelet(coeffR, coeffT, bookR);

% IMPORTANT:
% DLMAIM expects magnitudes (as in your original code)
dlmaimScore = dlmaim(abs(Rcoeff), abs(coeffR), abs(Acoeff), bookR);

% -------------------------
% HaarPSI (optional, for further improvement of the metric)
% -------------------------
%{
if exist('HaarPSI', 'file') == 2
    try
        haarScore = HaarPSI(refY, testY);
    catch ME
        haarScore = 1;
        warning('HaarPSI failed due to internal size inconsistency in the public implementation: %s', ME.message);
    end
end
%}
% -------------------------
% Final score
% -------------------------
%score = simProduct * dlmaimScore * haarScore;
score = simProduct * dlmaimScore;
end

% -------------------------
% Helper
% -------------------------
function out = rgb2gray_uint8safe(img)
if isinteger(img)
    img = double(img);
end
out = 0.299 * img(:,:,1) + 0.587 * img(:,:,2) + 0.114 * img(:,:,3);
end
