function results = TID2013YIQ(dataPath)
% TID2013YIQ  Evaluate dataset (TID2013-like) using the main quality function.
%   results = TID2013YIQ(dataPath)
%
% Input:
%   dataPath - root folder where dataset is stored.
%              Expected structure (default if missing):
%                dataPath/reference_images/I01.BMP ... I25.BMP
%                dataPath/distorted_images/Ixx_yy_z.BMP  (your prior folder convention)
%   If dataPath is not provided, current folder is used and the user will be asked
%   to set valid paths inside the script.
%
% Output:
%   results - structure with fields:
%     TidScores  - vector of all computed scores
%     mld        - matrix of luminance differences per distortion
%     CorADM     - per-distortion Spearman correlations (1x24)
%
% Notes:
%   - This script does not bundle TID2013 dataset. Provide the dataset locally
%     and set dataPath accordingly.
%   - This script accepts the CCS detection logic that was in original code:
%     for moderate luminance change (0.9 < LumDiff < 1.4) it converts to Qcom
%     chroma-like component before feeding to main().

if nargin < 1 || isempty(dataPath)
    dataPath = uigetdir(pwd, 'C:\Users\mobinam\Downloads\tid2013');
    if isequal(dataPath, 0)
        error('TID2013YIQ: no data path provided. Cancelled by user.');
    end
end

% Expected subfolders (adjust if your dataset is organized differently)
refFolder = fullfile(dataPath, 'reference_images');
distFolder = fullfile(dataPath, 'distorted_images');
mosFile = fullfile(dataPath, 'mos.txt');

% Validate existence
if ~isfolder(refFolder) || ~isfolder(distFolder)
    error('TID2013YIQ: expected folders not found. Make sure %s and %s exist.', refFolder, distFolder);
end
if ~isfile(mosFile)
    warning('TID2013YIQ: mos.txt not found at %s. Continue without MOS loading.', mosFile);
end

% Preallocate (TID2013 typically 25 refs * 24 distortions * 5 levels = 3000)
Tid2013 = zeros(3000, 1);
mld = zeros(24, 25*5);

iPoint = 0;
for iRef = 1:25
    % Reference filename pattern I01.BMP ... I25.BMP
    imNameRef = sprintf('I%02d.BMP', iRef);
    refPath = fullfile(refFolder, imNameRef);
    if ~isfile(refPath)
        error('Reference file missing: %s', refPath);
    end
    s = double(imread(refPath));
    lumR = 0.299 * s(:,:,1) + 0.587 * s(:,:,2) + 0.114 * s(:,:,3);

    for iDis = 1:24
        for iLevel = 1:5
            imNameDis = sprintf('I%02d_%02d_%d.BMP', iRef, iDis, iLevel); % adjust pattern if needed
            tPath = fullfile(distFolder, imNameDis);
            if ~isfile(tPath)
                % Try fallback naming convention used previously:
                tPath = fullfile(distFolder, sprintf('I%02d_%02d_%d.BMP', iRef, iDis, iLevel));
                if ~isfile(tPath)
                    error('Distorted file missing: %s (or other expected naming)', tPath);
                end
            end
            t = double(imread(tPath));
            lumT = 0.299 * t(:,:,1) + 0.587 * t(:,:,2) + 0.114 * t(:,:,3);

            iPoint = iPoint + 1;
            LumDiff = max(lumT(:) - lumR(:));
            mld(iDis, (iRef-1)*5 + iLevel) = LumDiff;

            % If moderate luminance change, compute Q-like chroma and use that
            if (0.9 < LumDiff) && (LumDiff < 1.4)
                QcomR = 0.211 * s(:,:,1) - 0.523 * s(:,:,2) + 0.312 * s(:,:,3);
                QcomT = 0.211 * t(:,:,1) - 0.523 * t(:,:,2) + 0.312 * t(:,:,3);
                Tid2013(iPoint) = main(QcomR, QcomT);
            else
                Tid2013(iPoint) = main(lumR, lumT);
            end
        end
    end
end

results.TidScores = Tid2013
results.mld = mld;

% Try to load MOS / subjective scores if present
if isfile(mosFile)
    try
        SB = load(mosFile); % expects plain text vector; adapt if format different
        results.MOS = SB;
        disp('MOS loaded.');
        % compute correlation
        spCorr = corr(Tid2013, SB, 'type', 'Spearman');
        fprintf('Overall Spearman correlation: %.4f\n', spCorr);
        % Also compute per-distortion correlations
        ADMGrDis = reshape(Tid2013, [5*25, 24]);
        mosMat = reshape(SB, [5*25, 24]);
        for j = 1:24
            CorADM(j) = corr(mosMat(:,j), ADMGrDis(:,j), 'type', 'Spearman');
        end
        results.CorADM = CorADM;
        figure; bar(CorADM); ylim([0 1]); xlabel('Distortion'); ylabel('Spearman');
    catch ME
        warning('Could not load/parse mos.txt: %s', ME.message);
    end
else
    warning('MOS file not found. Only raw scores saved in results.TidScores.');
end

% Save results to current folder
save('Tid2013_results.mat', 'results', '-v7.3');

end
