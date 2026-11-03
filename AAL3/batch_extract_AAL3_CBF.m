
%% =========================================================
% Batch AAL3 Gray Matter CBF Extraction
%
% Software: MATLAB + SPM25
% Study: eTRE ASL
%
% Input:
%   1. AAL3_batch_summary.csv
%   2. Individual meanCBF_T1coreg.nii
%   3. Individual AAL3_GM_CBFspace.nii
%
% IMPORTANT:
%   - No original AAL3 atlas required
%   - No spatial registration performed
%   - No SPM jobman initialization
%   - No minimum ROI voxel threshold
%   - Missing ROI = NaN
%   - Zero/negative finite CBF values retained and reported
%
% Output:
%   AAL3_CBF_mean_wide.csv
%   AAL3_CBF_valid_voxels_wide.csv
%   AAL3_CBF_GM_voxels_wide.csv
%   AAL3_CBF_long.csv
%   AAL3_CBF_extraction_QC.csv
%   AAL3_ROI_ID_columns.csv
%
%% =========================================================

clearvars;
clc;

%% =========================================================
% 1. Configure MATLAB and SPM25
%% =========================================================

workDir = '/Users/mac/AAL3';
spmDir = '/Users/mac/spm';

assert(isfolder(workDir), ...
    'AAL3 working directory not found.');

assert(isfile(fullfile(spmDir,'spm.m')), ...
    'Independent SPM installation not found.');

cd(workDir);

% Remove ExploreASL modified SPM path
exploreSPM = ...
    '/Users/mac/MatlabTools/ExploreASL-1.11.0/External/SPMmodified';

if contains(path,exploreSPM)
    rmpath(genpath(exploreSPM));
end

% Add independent SPM25
addpath(spmDir,'-begin');
rehash;

assert(strcmp(spm('Ver'),'SPM25'), ...
    'Current SPM is not SPM25.');

fprintf('\nSPM version: %s\n',spm('Ver'));
fprintf('SPM directory: %s\n',spm('Dir'));

% NOTE:
% Do not call spm_jobman('initcfg').
% This extraction requires only spm_vol and spm_read_vols.

%% =========================================================
% 2. Define study directories
%% =========================================================

aslDir = '/Users/mac/EXP/eTRE/ASL/asl_run';

summaryFile = fullfile(aslDir, ...
    'AAL3_batch_summary.csv');

assert(isfolder(aslDir), ...
    'ASL study directory not found.');

assert(isfile(summaryFile), ...
    'AAL3 batch summary not found.');

% Separate output directory
outputDir = fullfile(aslDir,'AAL3_CBF_Results');

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% =========================================================
% 3. Read successfully processed subjects
%% =========================================================

% Preserve subject IDs as strings
opts = detectImportOptions(summaryFile, ...
    'TextType','string');

opts = setvartype(opts,'Subject','string');
opts = setvartype(opts,'Status','string');

batchTable = readtable(summaryFile,opts);

status = upper(strtrim(string(batchTable.Status)));

successMask = status == "SUCCESS";

subjects = strtrim( ...
    string(batchTable.Subject(successMask)));

% Remove duplicate subjects, retaining original order
subjects = unique(subjects,'stable');

nSubjects = numel(subjects);

assert(nSubjects > 0, ...
    'No SUCCESS subjects found.');

fprintf('\n========================================\n');
fprintf('AAL3 GM-CBF EXTRACTION\n');
fprintf('========================================\n');

fprintf('Subjects selected: %d\n',nSubjects);

%% =========================================================
% 4. Define fixed AAL3 ROI labels
%% =========================================================

% Original AAL3 label range: 1-170
%
% Unused IDs:
% 35, 36, 81, 82
%
% Total valid anatomical labels: 166
%
% These label IDs only define output columns.
% They are NOT used for spatial registration.

roiLabels = setdiff( ...
    (1:170)', ...
    [35;36;81;82]);

nROI = numel(roiLabels);

assert(nROI == 166, ...
    'Unexpected number of AAL3 ROIs.');

% Fixed ROI column names
roiNames = cellstr(compose('ROI_%03d',roiLabels));

fprintf('Fixed AAL3 ROI labels: %d\n',nROI);

%% =========================================================
% 5. Initialize result matrices
%% =========================================================

% Mean CBF
MeanCBF = nan(nSubjects,nROI);

% Median CBF
MedianCBF = nan(nSubjects,nROI);

% Standard deviation
SDCBF = nan(nSubjects,nROI);

% Total GM ROI voxels
GMvoxels = zeros(nSubjects,nROI);

% Voxels with finite CBF values
ValidVoxels = zeros(nSubjects,nROI);

% Voxels with non-finite CBF values
InvalidVoxels = zeros(nSubjects,nROI);

% Finite zero CBF voxels
ZeroVoxels = zeros(nSubjects,nROI);

% Finite negative CBF voxels
NegativeVoxels = zeros(nSubjects,nROI);

%% Subject-level QC

qcStatus = strings(nSubjects,1);
qcMessage = strings(nSubjects,1);

qcSurvivingROIs = zeros(nSubjects,1);
qcValidROIs = zeros(nSubjects,1);

qcGMvoxels = zeros(nSubjects,1);
qcValidVoxels = zeros(nSubjects,1);
qcInvalidVoxels = zeros(nSubjects,1);

qcZeroVoxels = zeros(nSubjects,1);
qcNegativeVoxels = zeros(nSubjects,1);

%% =========================================================
% 6. Batch extract CBF
%% =========================================================

for s = 1:nSubjects

    sub = char(subjects(s));

    fprintf('\n----------------------------------------\n');
    fprintf('Subject %d/%d: %s\n', ...
        s,nSubjects,sub);
    fprintf('----------------------------------------\n');

    try

        %% 6.1 Input files

        cbfDir = fullfile(aslDir,sub, ...
            'FUNC','CBF_T1coreg');

        cbfFile = fullfile(cbfDir, ...
            'meanCBF_T1coreg.nii');

        roiFile = fullfile(cbfDir, ...
            'AAL3','AAL3_GM_CBFspace.nii');

        assert(isfile(cbfFile), ...
            'CBF file not found.');

        assert(isfile(roiFile), ...
            'Individual AAL3 GM atlas not found.');

        %% 6.2 Read images

        Vcbf = spm_vol(cbfFile);
        Vroi = spm_vol(roiFile);

        assert(numel(Vcbf)==1 && numel(Vroi)==1, ...
            'CBF and ROI must be single-volume 3D images.');

        %% 6.3 Check spatial geometry

        assert(isequal(Vcbf.dim,Vroi.dim), ...
            'CBF and ROI dimensions do not match.');

        affineDifference = ...
            max(abs(Vcbf.mat(:)-Vroi.mat(:)));

        assert(affineDifference < 1e-4, ...
            'CBF and ROI affine matrices do not match.');

        %% 6.4 Load voxel data

        CBF = spm_read_vols(Vcbf);
        ROI = spm_read_vols(Vroi);

        %% 6.5 Validate ROI labels

        assert(all(isfinite(ROI(:))), ...
            'ROI atlas contains NaN or Inf.');

        assert(all(abs(ROI(:)-round(ROI(:))) < 1e-5), ...
            'ROI atlas contains non-integer labels.');

        assert(all(ROI(:)>=0 & ROI(:)<=170), ...
            'ROI label outside expected range.');

        observedLabels = unique(ROI(ROI>0));

        unexpected = setdiff(observedLabels,roiLabels);

        assert(isempty(unexpected), ...
            'ROI atlas contains unexpected labels.');

        %% 6.6 Identify valid CBF values

        finiteMask = isfinite(CBF);

        %% 6.7 Extract every anatomical ROI

        for k = 1:nROI

            currentROI = roiLabels(k);

            % Individual GM-restricted ROI
            idx = (ROI == currentROI);

            % Number of anatomical GM voxels
            nGM = nnz(idx);

            GMvoxels(s,k) = nGM;

            % Missing ROI: leave CBF as NaN
            if nGM == 0
                continue;
            end

            % Extract CBF values in the ROI
            valuesAll = CBF(idx);

            % Finite CBF values only
            values = valuesAll(isfinite(valuesAll));

            nValid = numel(values);

            nInvalid = nGM-nValid;

            ValidVoxels(s,k) = nValid;
            InvalidVoxels(s,k) = nInvalid;

            % Missing valid CBF: leave NaN
            if nValid == 0
                continue;
            end

            %% Count special CBF values

            ZeroVoxels(s,k) = nnz(values==0);

            NegativeVoxels(s,k) = nnz(values<0);

            %% CBF statistics

            MeanCBF(s,k) = mean(values);

            MedianCBF(s,k) = median(values);

            % Standard deviation
            if nValid >= 2
                SDCBF(s,k) = std(values,0);
            else
                % SD not estimable from one voxel
                SDCBF(s,k) = NaN;
            end

        end

        %% 6.8 Subject quality control

        qcSurvivingROIs(s) = ...
            nnz(GMvoxels(s,:)>0);

        qcValidROIs(s) = ...
            nnz(ValidVoxels(s,:)>0);

        qcGMvoxels(s) = ...
            sum(GMvoxels(s,:));

        qcValidVoxels(s) = ...
            sum(ValidVoxels(s,:));

        qcInvalidVoxels(s) = ...
            sum(InvalidVoxels(s,:));

        qcZeroVoxels(s) = ...
            sum(ZeroVoxels(s,:));

        qcNegativeVoxels(s) = ...
            sum(NegativeVoxels(s,:));

        qcStatus(s) = "SUCCESS";
        qcMessage(s) = "Completed";

        fprintf('Surviving GM ROIs: %d\n', ...
            qcSurvivingROIs(s));

        fprintf('ROIs with valid CBF: %d\n', ...
            qcValidROIs(s));

        fprintf('GM ROI voxels: %d\n', ...
            qcGMvoxels(s));

        fprintf('Valid CBF voxels: %d\n', ...
            qcValidVoxels(s));

        fprintf('Zero CBF voxels: %d\n', ...
            qcZeroVoxels(s));

        fprintf('Negative CBF voxels: %d\n', ...
            qcNegativeVoxels(s));

        fprintf('SUCCESS: %s\n',sub);

    catch ME

        qcStatus(s) = "FAILED";
        qcMessage(s) = string(ME.message);

        % Reset failed subject data
        MeanCBF(s,:) = NaN;
        MedianCBF(s,:) = NaN;
        SDCBF(s,:) = NaN;

        GMvoxels(s,:) = 0;
        ValidVoxels(s,:) = 0;
        InvalidVoxels(s,:) = 0;
        ZeroVoxels(s,:) = 0;
        NegativeVoxels(s,:) = 0;

        fprintf(2,'FAILED: %s\n',sub);
        fprintf(2,'%s\n',ME.message);

    end

end

%% =========================================================
% 7. Identify successfully extracted subjects
%% =========================================================

finalSuccess = qcStatus=="SUCCESS";

successSubjects = subjects(finalSuccess);

fprintf('\n========================================\n');
fprintf('EXTRACTION SUMMARY\n');
fprintf('========================================\n');

fprintf('Selected subjects: %d\n',nSubjects);

fprintf('Successful: %d\n',nnz(finalSuccess));

fprintf('Failed: %d\n',nnz(~finalSuccess));

%% =========================================================
% 8. Wide-format mean CBF table
%% =========================================================

meanTable = array2table( ...
    MeanCBF(finalSuccess,:), ...
    'VariableNames',roiNames);

meanTable = addvars( ...
    meanTable, ...
    successSubjects, ...
    'Before',1, ...
    'NewVariableNames','Subject');

meanFile = fullfile(outputDir, ...
    'AAL3_CBF_mean_wide.csv');

writetable(meanTable,meanFile);

%% =========================================================
% 9. Wide-format valid voxel count table
%% =========================================================

validTable = array2table( ...
    ValidVoxels(finalSuccess,:), ...
    'VariableNames',roiNames);

validTable = addvars( ...
    validTable, ...
    successSubjects, ...
    'Before',1, ...
    'NewVariableNames','Subject');

validFile = fullfile(outputDir, ...
    'AAL3_CBF_valid_voxels_wide.csv');

writetable(validTable,validFile);

%% =========================================================
% 10. Wide-format GM voxel count table
%% =========================================================

gmTable = array2table( ...
    GMvoxels(finalSuccess,:), ...
    'VariableNames',roiNames);

gmTable = addvars( ...
    gmTable, ...
    successSubjects, ...
    'Before',1, ...
    'NewVariableNames','Subject');

gmFile = fullfile(outputDir, ...
    'AAL3_CBF_GM_voxels_wide.csv');

writetable(gmTable,gmFile);

%% =========================================================
% 11. Long-format CBF table
%% =========================================================

% Select successful subjects
Mmean = MeanCBF(finalSuccess,:);
Mmedian = MedianCBF(finalSuccess,:);
Msd = SDCBF(finalSuccess,:);

Mgm = GMvoxels(finalSuccess,:);
Mvalid = ValidVoxels(finalSuccess,:);
Minvalid = InvalidVoxels(finalSuccess,:);
Mzero = ZeroVoxels(finalSuccess,:);
Mnegative = NegativeVoxels(finalSuccess,:);

nSuccess = numel(successSubjects);

% Subject-major long format
longSubject = repelem(successSubjects,nROI);

longROI = repmat(roiLabels,nSuccess,1);

% Transpose before reshape to preserve Subject x ROI order
longMean = reshape(Mmean.',[],1);
longMedian = reshape(Mmedian.',[],1);
longSD = reshape(Msd.',[],1);

longGM = reshape(Mgm.',[],1);
longValid = reshape(Mvalid.',[],1);
longInvalid = reshape(Minvalid.',[],1);
longZero = reshape(Mzero.',[],1);
longNegative = reshape(Mnegative.',[],1);

longTable = table( ...
    longSubject, ...
    longROI, ...
    longGM, ...
    longValid, ...
    longInvalid, ...
    longZero, ...
    longNegative, ...
    longMean, ...
    longMedian, ...
    longSD, ...
    'VariableNames',{ ...
    'Subject', ...
    'ROI_ID', ...
    'GM_Voxels', ...
    'CBF_Valid_Voxels', ...
    'CBF_Invalid_Voxels', ...
    'CBF_Zero_Voxels', ...
    'CBF_Negative_Voxels', ...
    'Mean_CBF', ...
    'Median_CBF', ...
    'SD_CBF'});

longFile = fullfile(outputDir, ...
    'AAL3_CBF_long.csv');

writetable(longTable,longFile);

%% =========================================================
% 12. Subject-level extraction QC table
%% =========================================================

qcTable = table( ...
    subjects, ...
    qcStatus, ...
    qcMessage, ...
    qcSurvivingROIs, ...
    qcValidROIs, ...
    qcGMvoxels, ...
    qcValidVoxels, ...
    qcInvalidVoxels, ...
    qcZeroVoxels, ...
    qcNegativeVoxels, ...
    'VariableNames',{ ...
    'Subject', ...
    'Status', ...
    'Message', ...
    'Surviving_ROIs', ...
    'ROIs_with_valid_CBF', ...
    'GM_ROI_Voxels', ...
    'Valid_CBF_Voxels', ...
    'Invalid_CBF_Voxels', ...
    'Zero_CBF_Voxels', ...
    'Negative_CBF_Voxels'});

qcFile = fullfile(outputDir, ...
    'AAL3_CBF_extraction_QC.csv');

writetable(qcTable,qcFile);

%% =========================================================
% 13. ROI column mapping
%% =========================================================

roiMap = table( ...
    roiLabels(:), ...
    string(roiNames(:)), ...
    'VariableNames',{ ...
    'ROI_ID','Column_Name'});

roiMapFile = fullfile(outputDir, ...
    'AAL3_ROI_ID_columns.csv');

writetable(roiMap,roiMapFile);

%% =========================================================
% 14. Final QC summary
%% =========================================================

fprintf('\n========================================\n');
fprintf('FINAL AAL3 GM-CBF RESULTS\n');
fprintf('========================================\n');

fprintf('Successful subjects: %d\n',nSuccess);

fprintf('Fixed ROI columns: %d\n',nROI);

fprintf('Long-format rows: %d\n',height(longTable));

% Missing ROI-CBF entries
missingCBF = sum(isnan(Mmean(:)));

totalEntries = numel(Mmean);

fprintf('Missing ROI-CBF entries: %d / %d\n', ...
    missingCBF,totalEntries);

fprintf('Missing fraction: %.2f%%\n', ...
    100*missingCBF/totalEntries);

% ROI count statistics
fprintf('Mean surviving ROIs: %.2f\n', ...
    mean(qcSurvivingROIs(finalSuccess)));

fprintf('Minimum surviving ROIs: %d\n', ...
    min(qcSurvivingROIs(finalSuccess)));

fprintf('Maximum surviving ROIs: %d\n', ...
    max(qcSurvivingROIs(finalSuccess)));

%% =========================================================
% 15. Output file information
%% =========================================================

fprintf('\nOutput directory:\n%s\n',outputDir);

fprintf('\nGenerated files:\n');

fprintf('1. AAL3_CBF_mean_wide.csv\n');
fprintf('2. AAL3_CBF_valid_voxels_wide.csv\n');
fprintf('3. AAL3_CBF_GM_voxels_wide.csv\n');
fprintf('4. AAL3_CBF_long.csv\n');
fprintf('5. AAL3_CBF_extraction_QC.csv\n');
fprintf('6. AAL3_ROI_ID_columns.csv\n');

fprintf('\n========================================\n');
fprintf('AAL3 CBF EXTRACTION COMPLETED\n');
fprintf('========================================\n');
