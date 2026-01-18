function [R, A] = decouple_wavelet(O, T, S)
% DECOUPLE_WAVELET
% Wavelet-domain decoupling of detail losses and additive impairments (DPD).
%
% Inputs:
%   O - reference wavelet coefficient vector (from wavedec2)
%   T - test/distorted wavelet coefficient vector (from wavedec2)
%   S - bookkeeping matrix from wavedec2
%
% Outputs:
%   R - restored/detail-related coefficients
%   A - additive-impairment coefficients
%
% NOTE:
%   This implementation follows the ADM-style decoupling strategy used in
%   public image quality assessment code. This code is not originally authored
%   by the repository owner and is included for reproducibility.
%
%   The decoupling is done directly in the flattened coefficient vector using
%   explicit subband indexing. This avoids detcoef2/wrcoef2 incompatibilities.

% Initialize
R = zeros(size(O));

% -------------------------
% Copy approximation band
% -------------------------
% Approximation band corresponds to (lambda=4, theta=1)
[s, e] = getSubbandSE(4, 1, S);
R(s:e) = T(s:e);

% -------------------------
% Process detail subbands
% -------------------------
for lambda = 1:4
    for theta = 2:4   % H, V, D subbands

        [s, e] = getSubbandSE(lambda, theta, S);

        Osub = O(s:e);
        Tsub = T(s:e);

        % Gain factor (detail loss), clamped to [0, 1]
        k = Tsub ./ (Osub + 1e-30);
        k = min(max(k, 0), 1);

        Rsub = k .* Osub;

        % Special case: contrast enhancement correction
        if theta == 2 || theta == 3
            % Horizontal & vertical used for orientation check
            if theta == 2
                [sV, eV] = getSubbandSE(lambda, 3, S);
                OV = O(sV:eV);
                TV = T(sV:eV);
                OA = calAngle(Osub, OV);
                TA = calAngle(Tsub, TV);
            else
                [sH, eH] = getSubbandSE(lambda, 2, S);
                OH = O(sH:eH);
                TH = T(sH:eH);
                OA = calAngle(OH, Osub);
                TA = calAngle(TH, Tsub);
            end

            diff = abs(OA - TA) * 180 / pi;
            pos = diff < 1;

            Rsub(pos) = Tsub(pos);
        end

        R(s:e) = Rsub;
    end
end

% -------------------------
% Additive impairment
% -------------------------
A = T - R;

end

% ==========================================================
% Helper functions
% ==========================================================

function [s, e] = getSubbandSE(lambda, theta, S)
% GETSUBBANDSE
% Returns start/end indices of a subband in wavedec2 coefficient vector.

if lambda == 4 && theta == 1
    s = 1;
    e = S(1,1) * S(1,2);
    return;
end

lambda2 = 5 - lambda;
s = S(1,1) * S(1,2) + 1;

for i = 2:lambda2
    s = s + 3 * S(i,1) * S(i,2);
end

theta2 = theta - 2;
s = s + theta2 * S(lambda2+1,1) * S(lambda2+1,2);
e = s + S(lambda2+1,1) * S(lambda2+1,2) - 1;
end

function A = calAngle(H, V)
% CALANGLE Compute orientation angle from horizontal and vertical components
A = atan(V ./ (H + 1e-30)) + pi * (H < 0);
end

