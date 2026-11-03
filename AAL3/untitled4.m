
%% =========================================================
% AAL3 124 ROI CBF - Final Excel Export
% Subjects: 130
% Bilateral pairs: 59
% Midline ROIs: 6
% Individual voxel threshold: >= 10
%% =========================================================

clear;
clc;

%% 1. Paths

aslDir = '/Users/mac/EXP/eTRE/ASL/asl_run';

dataDir = fullfile(aslDir,'AAL3_CBF_Results');

roiFile = '/Users/mac/AAL3/AAL3_ROI_124_retained.txt';

meanFile = fullfile(dataDir,'AAL3_CBF_mean_wide.csv');
validFile = fullfile(dataDir,'AAL3_CBF_valid_voxels_wide.csv');
gmFile = fullfile(dataDir,'AAL3_CBF_GM_voxels_wide.csv');

excelFile = fullfile(dataDir,'AAL3_CBF_124ROI_Final.xlsx');

minVoxels = 10;

%% 2. Check files

assert(isfile(roiFile),'124 ROI dictionary not found.');
assert(isfile(meanFile),'Mean CBF CSV not found.');
assert(isfile(validFile),'Valid voxel CSV not found.');
assert(isfile(gmFile),'GM voxel CSV not found.');

%% 3. Read 124 ROI names

lines = readlines(roiFile);

roiIDs = [];
roiNames = strings(0,1);

for i = 1:numel(lines)

    t = regexp(strtrim(char(lines(i))), ...
        '^(\d+)\s+([A-Za-z0-9_]+)$','tokens','once');

    if ~isempty(t)
        roiIDs(end+1,1) = str2double(t{1});
        roiNames(end+1,1) = string(t{2});
    end
end

assert(numel(roiIDs)==124,'ROI count must equal 124.');
assert(numel(unique(roiIDs))==124,'Duplicate ROI IDs.');

nROI = numel(roiIDs);

roiCols = cellstr(compose('ROI_%03d',roiIDs));

%% 4. Read original CBF data

Tmean = readtable(meanFile,'TextType','string');
Tvalid = readtable(validFile,'TextType','string');
Tgm = readtable(gmFile,'TextType','string');

subjects = strtrim(string(Tmean.Subject));
nSubjects = numel(subjects);

[ok1,idx1] = ismember(subjects,strtrim(string(Tvalid.Subject)));
[ok2,idx2] = ismember(subjects,strtrim(string(Tgm.Subject)));

assert(all(ok1) && all(ok2),'Subject matching failed.');

assert(numel(unique(subjects))==nSubjects, ...
    'Duplicate subject IDs.');

% Check all required ROI columns

assert(all(ismember(roiCols,Tmean.Properties.VariableNames)));
assert(all(ismember(roiCols,Tvalid.Properties.VariableNames)));
assert(all(ismember(roiCols,Tgm.Properties.VariableNames)));

MeanRaw = table2array(Tmean(:,roiCols));
ValidVoxels = table2array(Tvalid(idx1,roiCols));
GMVoxels = table2array(Tgm(idx2,roiCols));

assert(all(isfinite(ValidVoxels(:))) && ...
    all(ValidVoxels(:)>=0) && ...
    all(ValidVoxels(:)==round(ValidVoxels(:))), ...
    'Invalid voxel count data.');

%% 5. Identify bilateral pairs

pairL = [];
pairR = [];
paired = false(nROI,1);

for k = 1:nROI

    token = regexp(char(roiNames(k)), ...
        '^(.*)_L$','tokens','once');

    if isempty(token)
        continue;
    end

    rightName = string(token{1}) + "_R";

    j = find(roiNames==rightName,1);

    if ~isempty(j)
        pairL(end+1,1) = k;
        pairR(end+1,1) = j;
        paired([k j]) = true;
    end
end

midline = find(~paired);

nPairs = numel(pairL);
nMidline = numel(midline);

fprintf('\nBilateral pairs: %d\n',nPairs);
fprintf('Independent ROIs: %d\n',nMidline);

% Verify all 124 ROIs are correctly classified

assert(nPairs==59 && nMidline==6 && ...
    2*nPairs+nMidline==124, ...
    'ROI pairing error. Check the ROI dictionary.');

% Check that every lateral ROI has a counterpart

for k = 1:nROI

    name = char(roiNames(k));

    if endsWith(name,'_L') || endsWith(name,'_R')
        assert(paired(k), ...
            'Unpaired lateral ROI: %s',name);
    end
end

%% 6. Apply individual-level voxel threshold

MeanFinal = MeanRaw;

pairFailed = false(nSubjects,nPairs);

% Bilateral ROIs: both sides must pass

for p = 1:nPairs

    L = pairL(p);
    R = pairR(p);

    bad = ValidVoxels(:,L)<minVoxels | ...
          ValidVoxels(:,R)<minVoxels;

    MeanFinal(bad,[L R]) = NaN;

    pairFailed(:,p) = bad;
end

% Midline ROIs: assessed independently

for k = midline'

    bad = ValidVoxels(:,k)<minVoxels;

    MeanFinal(bad,k) = NaN;
end

%% 7. Generate output tables

Tfinal = array2table(MeanFinal,'VariableNames',roiCols);
Tfinal = addvars(Tfinal,subjects,'Before',1, ...
    'NewVariableNames','Subject');

Traw = array2table(MeanRaw,'VariableNames',roiCols);
Traw = addvars(Traw,subjects,'Before',1, ...
    'NewVariableNames','Subject');

Tvalid124 = array2table(ValidVoxels,'VariableNames',roiCols);
Tvalid124 = addvars(Tvalid124,subjects,'Before',1, ...
    'NewVariableNames','Subject');

Tgm124 = array2table(GMVoxels,'VariableNames',roiCols);
Tgm124 = addvars(Tgm124,subjects,'Before',1, ...
    'NewVariableNames','Subject');

Tdictionary = table(roiIDs,roiNames, ...
    string(roiCols(:)), ...
    'VariableNames',{'ROI_ID','ROI_Name','Column_Name'});

Tpairing = table( ...
    roiIDs(pairL),roiNames(pairL), ...
    roiIDs(pairR),roiNames(pairR), ...
    'VariableNames',{'Left_ID','Left_Name','Right_ID','Right_Name'});

Tqc = table( ...
    subjects, ...
    sum(isfinite(MeanFinal),2), ...
    sum(~isfinite(MeanFinal),2), ...
    sum(pairFailed,2), ...
    'VariableNames',{'Subject','Valid_ROIs', ...
    'Missing_ROIs','Failed_Bilateral_Pairs'});

% Cohort-level ROI quality control

Tcohort = table( ...
    roiIDs,roiNames, ...
    median(ValidVoxels,1)', ...
    mean(ValidVoxels,1)', ...
    sum(ValidVoxels>0,1)', ...
    'VariableNames',{'ROI_ID','ROI_Name', ...
    'Median_Voxels','Mean_Voxels','Subjects_With_Data'});

%% 8. Validate bilateral missingness

for p = 1:nPairs

    assert(isequal( ...
        isnan(MeanFinal(:,pairL(p))), ...
        isnan(MeanFinal(:,pairR(p)))), ...
        'Bilateral missingness mismatch.');
end

%% 9. Export Excel safely

% Write to a temporary workbook first.
% Replace the previous final workbook only after successful export.

tempFile = fullfile(dataDir,'AAL3_CBF_124ROI_temp.xlsx');

if isfile(tempFile)
    delete(tempFile);
end

writetable(Tfinal,tempFile,'Sheet','Mean_CBF_Final');
writetable(Traw,tempFile,'Sheet','Mean_CBF_Raw');
writetable(Tvalid124,tempFile,'Sheet','Valid_Voxels');
writetable(Tgm124,tempFile,'Sheet','GM_Voxels');
writetable(Tdictionary,tempFile,'Sheet','ROI_Dictionary');
writetable(Tpairing,tempFile,'Sheet','Bilateral_Pairs');
writetable(Tqc,tempFile,'Sheet','Subject_QC');
writetable(Tcohort,tempFile,'Sheet','Cohort_ROI_QC');

assert(isfile(tempFile),'Excel export failed.');

% Replace previous 124 ROI Excel
if isfile(excelFile)
    try
        delete(excelFile);
    catch
        error('Please close the existing Excel workbook and rerun.');
    end
end

[success,msg] = movefile(tempFile,excelFile);

assert(success,'Cannot replace final Excel: %s',msg);

%% 10. Final report

fprintf('\n====================================\n');
fprintf('AAL3 CBF EXPORT COMPLETED\n');
fprintf('====================================\n');

fprintf('Subjects: %d\n',nSubjects);
fprintf('Retained ROIs: %d\n',nROI);
fprintf('Bilateral pairs: %d\n',nPairs);
fprintf('Midline ROIs: %d\n',nMidline);

fprintf('Mean valid ROIs/subject: %.2f\n', ...
    mean(sum(isfinite(MeanFinal),2)));

fprintf('Minimum valid ROIs: %d\n', ...
    min(sum(isfinite(MeanFinal),2)));

fprintf('Maximum valid ROIs: %d\n', ...
    max(sum(isfinite(MeanFinal),2)));

fprintf('Total missing values: %d\n', ...
    sum(isnan(MeanFinal(:))));

fprintf('Failed bilateral pairs: %d\n', ...
    sum(pairFailed(:)));

fprintf('\nExcel saved:\n%s\n',excelFile);

fprintf('\nOriginal 166 ROI CSV files unchanged.\n');

fprintf('\nALL CHECKS PASSED.\n');
