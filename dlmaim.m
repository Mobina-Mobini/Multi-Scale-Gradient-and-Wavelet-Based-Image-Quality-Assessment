function q = dlmaim(R3, O3, A3, S)
% DLMAIM
% Detail Loss Measure with Additive Impairment Modulation.
%
% Inputs:
%   R3 - absolute value of restored/detail coefficients
%   O3 - absolute value of reference coefficients
%   A3 - absolute value of additive impairment coefficients
%   S  - bookkeeping matrix from wavedec2
%
% Output:
%   q  - scalar quality score
%
% NOTE:
%   This implementation follows the ADM-style public code used in the
%   literature (e.g., Li et al., IEEE). It operates directly on the
%   wavedec2 coefficient vector using explicit subband indexing.

% -------------------------
% Part 1: Detail restoration term
% -------------------------
nume = 0;
deno = 0;

for lambda = 1:4
    for theta = 2:4   % detail subbands only

        [s, e] = getSubbandSE(lambda, theta, S);

        Rsub = reshape(R3(s:e), S(6-lambda,1), S(6-lambda,2));
        Osub = reshape(O3(s:e), S(6-lambda,1), S(6-lambda,2));

        % Remove border region (center crop)
        Rsub = centerCrop(Rsub);
        Osub = centerCrop(Osub);

        a = sum(Rsub(:).^3)^(1/3);
        b = sum(Osub(:).^3)^(1/3);

        nume = nume + a;
        deno = deno + b;
    end
end

q1 = nume / (deno + eps);

% -------------------------
% Part 2: Additive impairment penalty
% -------------------------
nume = 0;
deno = numel(A3);

for lambda = 1:4
    for theta = 2:4

        [s, e] = getSubbandSE(lambda, theta, S);

        Asub = reshape(A3(s:e), S(6-lambda,1), S(6-lambda,2));
        Asub = centerCrop(Asub);

        SD = std(Asub(:));
        a = sum((Asub(:) + SD).^3)^(1/3);

        nume = nume + a;
    end
end

a1 = 10;   % tuning constant (kept from original code)
q2 = nume / (deno + eps);

% -------------------------
% Final DLMAIM score
% -------------------------
q = q1 - a1 * q2;

end

% ==========================================================
% Helper functions
% ==========================================================

function B = centerCrop(A)
% CENTERCROP Remove border region to avoid boundary artifacts
marg = 0.1;  % 10% margin
[h, w] = size(A);
mh = max(1, round(h * marg));
mw = max(1, round(w * marg));
B = A(mh:h-mh+1, mw:w-mw+1);
end

function [s, e] = getSubbandSE(lambda, theta, S)
% GETSUBBANDSE Return start/end indices of a subband
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
