% This script converst the cato output files per subject (.mat) into a 
% single .mat object including regions descriptions.

%% set everything up
projPath = '/localpath';
rPath = '/usr/local/bin/Rscript';
reconMeth = 'csd_dti';
atlas = 'lausanne120'; % aparc lausanne120 schaefer100-yeo7
demFile = 'ITQ_demographics.csv';
hThr = 'threshold10';
nosIndex= 1 ;       % streamline matrix index
faIndex = 3 ;       % fa matrix index
mdIndex = 6;        % md matrix index
svdIndex = 13;
qThr = 3 ;          % number of IQR to define outliers
pThr = 0.7 ;        % percentage subjects that should have connection

figPath = [projPath '/figures/' atlas '/' hThr];
catoDir = [projPath '/data/' atlas,'_',reconMeth];

roiFile = [projPath '/data/' atlas '_regions_' hThr '.txt'];
lfsFile = [projPath '/data/',atlas '_thickness_lh.csv'];
matPattern = ['*connectivity_' reconMeth '_' atlas '.mat'];

addpath(genpath(projPath));
addpath('/local_path/Simple-Brain-Plot');
addpath('/local_path/cbrewer');


%% connectivity matrices, demographics and volume
matName = [projPath '/data/sandc_' reconMeth '_' atlas '.mat'];

if isfile(matName)
    fprintf('Loading %s \n', matName);
    load(matName, 'sandcMat');
else
    % load data
    sandcMat = mergeCATO(catoDir, matPattern);       % merge CATO matrices
    tableDemo = readtable([projPath '/data/' demFile]); % load demographics
    lhFS = readtable(lfsFile); % only need brain volume
    
    % combine into single object
    tableDemo = tableDemo(:, {'subject', 'insomnia_subtype', 'Age', 'Sex', 'study', 'ISI_score'});
    tableDemo.Properties.VariableNames = ...
        {'subject' 'subtype' 'age' 'male' 'study' 'isi'}; % names
    [L,i] = ismember(sandcMat.subjects, tableDemo.subject);  % match index
    sandcMat.subjects = sandcMat.subjects(L);          % subset
    sandcMat.connectivity = sandcMat.connectivity(:,:,:,L); % subset  
    sandcMat.table = tableDemo(i(i~=0),:);              % add ordered table
    
    Sub = regexp(lhFS{:,1},'sub-[a-z0-9]{8}', 'match'); % subject
    lhFS.subject = vertcat(Sub{:});
    [L,i] = ismember(lhFS.subject, sandcMat.subjects);  % match index
    lhFS = lhFS(L,:);
    sandcMat.table.BrainSegVolNotVent = lhFS{i(i~=0),'BrainSegVolNotVent'};
    
    % only store weights of interest to save memory
    sandcMat.connectivity = sandcMat.connectivity(:,:,[nosIndex, faIndex, mdIndex, svdIndex],:);
    sandcMat.weightDescriptions = sandcMat.weightDescriptions([nosIndex, faIndex, mdIndex, svdIndex]);
    
    save(matName, 'sandcMat');
    clear L i tableDemo matName Sub
end


%% quality check
% remove outliers using IQR
%   colum1:     mean streamlines, 
%   column2:    mean FA 
%   column3:    mean MD
%   column4:    average prevalence of present connections
%   column5:    average prevalence of connections not found in this subject 
%   column6:    number of connected regions

% new indexes
nosIndex= 1 ;       
faIndex = 2 ;      
mdIndex = 3;        
svdIndex = 4;

Mat = sandcMat.connectivity;
Mat(Mat==0) = NaN;
S = zeros(size(Mat, 4), 6); %pre-assign matrix

S(:,1) = nanmean(Mat(:,:,nosIndex,:), 1:2); 
S(:,2) = nanmean(Mat(:,:,faIndex,:), 1:2);
S(:,3) = nanmean(Mat(:,:,mdIndex,:), 1:2);

B = squeeze(Mat(:,:,nosIndex,:) > 0) ;          % binary matrix per subject
P = mean(Mat(:, :, nosIndex, :) > 0, 4);        % prevelance of each connection
S(:,4) = mean(B .* P, 1:2);                     % prevalance of present connections

C = B + eye(size(B(:,:,1)));                    % set B midline to 1
S(:,5) = mean(~C .* P, 1:2);                    % invert C mask for prevalance of abscent connections

S(:,6) = sum(Mat(:,:,nosIndex,:) > 0, 1:2) / 2 ;


% IQR
low = quantile(S, 0.25) - (qThr .* iqr(S)); % lower boundary
upp = quantile(S, 0.75) + (qThr .* iqr(S)); % upper boundary
outliers = S > upp | S < low ;              % outlier matrix
outliers = any(outliers, 2) ;    
    % strict: excluded if outlier in any measure


% remove outliers
cleanMat = sandcMat; % copy
cleanMat.subjects = cleanMat.subjects(~outliers);
cleanMat.connectivity = cleanMat.connectivity(:,:,:,~outliers);
cleanMat.table = cleanMat.table(~outliers,:);

% save outliers
writetable(sandcMat.table(outliers,:),[projPath, '/data/', atlas,'_outliers_', reconMeth,'.txt'],...
    'Delimiter', ' ')  

clear sandcMat Mat outliers P B C low upp


%% Table
cleanMat.table.subtype = categorical(cleanMat.table.subtype);
cleanMat.table.study = categorical(cleanMat.table.study);
cleanMat.table.male = categorical(cleanMat.table.male);


%% threshold using prevalence
Mat = cleanMat.connectivity;
SL = squeeze(Mat(:,:,faIndex,:));
P = mean(SL > 0, 3);
mask = double(P >= pThr); 
thrMat = Mat .* mask ;    
    clear P SL

    
%% Combat harmonization of FA and MD    
% ~ subtype + age + sex + brainvolume
% ercp is reference batch (best sequence and most subjects)
% harmonize across all edges in the brain

projT = cleanMat.table(:, {'subtype', 'age', 'male', 'BrainSegVolNotVent', 'study'});
  
% FA matrices
faMat = squeeze(thrMat(:,:,faIndex,:));
faCol = squareform3d(faMat); % transform (cols=edges and rows=subjects)
faCombat = neuroCombat_wrap(faCol, projT, rPath, [projPath, '/data']);
faMatCombat = squareform3d(faCombat); 

% MD matrices
mdMat = squeeze(thrMat(:,:,mdIndex,:));
mdCol = squareform3d(mdMat); % transform (cols=edges and rows=subjects)
mdCombat = neuroCombat_wrap(mdCol, projT, rPath, [projPath, '/data']);
mdMatCombat = squareform3d(mdCombat); 

% SVD matrices (streamline volume density)
svdMat = squeeze(thrMat(:,:,svdIndex,:));
svdCol = squareform3d(svdMat); % transform (cols=edges and rows=subjects)
svdCombat = neuroCombat_wrap(svdCol, projT, rPath, [projPath, '/data']);
svdMatCombat = squareform3d(svdCombat); 
    clearvars faCombat mdCombat svdCombat faCol mdCol svdCol

    
%% subset using insomnia-linked ROIs
clearvars Mat
Mat(:,:,1,:) = faMatCombat;
Mat(:,:,2,:) = mdMatCombat;
Mat(:,:,3,:) = svdMatCombat;

% Store whole brain combat data
cleanMat.combat.connectivity = Mat;
cleanMat.combat.regionDescriptions = cleanMat.regionDescriptions;
cleanMat.combat.weightDescriptions = [{'faCombat'};{'mdCombat'};{'svdCombat'}];

% copy of unharmonized data
Mat(:,:,4,:) = faMat;
Mat(:,:,5,:) = mdMat;
Mat(:,:,6,:) = svdMat;

% load regions of interest
roiNames = importdata(roiFile);
iROI = ismember(cleanMat.regionDescriptions, roiNames);

% make subset of whole brain harmonized and unharmonized data
cleanMat.NKE.weightDescriptions = [{'faCombat'};{'mdCombat'};{'svdCombat'}; ...
       cleanMat.weightDescriptions([faIndex, mdIndex,svdIndex])];
cleanMat.NKE.connectivity = Mat(iROI,iROI,:,:);

cleanMat.NKE.regionDescriptions = roiNames;
cleanMat.NKE.mask = mask(iROI,iROI);
cleanMat.NKE.threshold = pThr;
cleanMat.NKE.table = projT;
save([projPath '/data/' atlas '_' hThr '_' reconMeth '_cleanMat.mat'], 'cleanMat');