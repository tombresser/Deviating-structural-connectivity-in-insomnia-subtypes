clear

%% set everything up
rng(92) % reproducable seed

projPath = '/local_path';
reconMeth = 'csd_dti';
atlas = 'lausanne120';  % aparc lausanne120 schaefer100-yeo7
demFile = 'ITQ_demographics.csv';
hThr = 'threshold10';
nPerm = 10000;
combat = 'true';
brainvol = 'true';

result_path = [projPath '/figures/', atlas];
path_yeo = ['/local_path/parcellation/Yeo7_in_', atlas, '.mat'];

% plot settings
titles = ["Hetrogennious, Inssomnia Disorder", "highly distressed", ...
        "moderately distressed, reward sensitive", ...
        "moderately distressed, reward insensitive", ...
        "slightly distressed, high reactive", ...
        "slightly distressed, low reactive"];

% create path (names)
cachePath = [projPath '/data'];
if ~exist(result_path, 'dir')
    mkdir(result_path);
end

% prefix
prefix = [reconMeth,'_'];
if strcmp(combat, 'true')
    prefix = append('combat_', prefix);
end
if strcmp(brainvol, 'true')
    prefix = append('brainvol_', prefix);
end


    
% add paths
addpath(genpath(projPath));
addpath('/local_path/circularGraph');
addpath('local_path/cbrewer');
addpath('/local_path/Simple-Brain-Plot');
addpath('/local_path/subplot_tight');
addpath('/local_path/TolColors');

% load full matrices and subnetwork of interest
load([cachePath '/' atlas '_' hThr '_' reconMeth '_cleanMat.mat'], 'cleanMat');


%% Load cortical thickness and surface area results
% add date variable
names=["ID","ITQ1","ITQ2","ITQ3","ITQ4","ITQ5"];



%% prepare covariate matrices
% table
T = cleanMat.NKE.table;
T.subtype = double(T.subtype);
T.male = double(T.male);
T.study = double(T.study);

% zscore age and brainvolume
T.age = zscore(T.age);
T.brainVol = zscore(T.BrainSegVolNotVent);

% turn categorical into dummy
xMat = T(:,{'subtype','age','male' ,'brainVol', 'study'});
xMat = table2array(xMat);
d1=dummyvar(xMat(:,1));
d2=dummyvar(xMat(:,3));
d3=dummyvar(xMat(:,5));

if strcmp(brainvol, 'false')
    % variable matrices
    xITQ = [d1(:,2:end), xMat(:,2), d2(:,2)];
    xID = [any(d1(:,2:end),2), xMat(:,2), d2(:,2)];
elseif strcmp(brainvol, 'true')
    % include brainvol as covar
    xITQ = [d1(:,2:end), xMat(:,2), d2(:,2), xMat(:,4)];
    xID = [any(d1(:,2:end),2), xMat(:,2), d2(:,2), xMat(:,4)];
end

if strcmp(combat, 'false')
   xITQ(:,[end+1:end+3]) = d3(:,2:end);
   xID(:,[end+1:end+3]) = d3(:,2:end);
end

clearvars xMat d1 d2


%% prepare FA and MD matrices

% select correct connectome matrices
if strcmp(combat, 'true')
    FA = squeeze(cleanMat.NKE.connectivity(:,:,1,:));
    MD = squeeze(cleanMat.NKE.connectivity(:,:,2,:));
    SVD = squeeze(cleanMat.NKE.connectivity(:,:,3,:)); 
elseif strcmp(combat, 'false')
    FA = squeeze(cleanMat.NKE.connectivity(:,:,4,:));
    MD = squeeze(cleanMat.NKE.connectivity(:,:,5,:));
    SVD = squeeze(cleanMat.NKE.connectivity(:,:,6,:)); 
end

% yMat with left side (every column is an edge)
FA(isnan(FA))=0; % ugly solution
FA_col = squareform3d(FA); % edge2column
FA_col(FA_col==0)=nan;     % ugly solution

% yMat with left side (every column is an edge)
MD(isnan(MD))=0;           % ugly solution
MD_col = squareform3d(MD); % edge2column
MD_col(MD_col==0)=nan;     % ugly solution

% yMat with left side (every column is an edge)
SVD(isnan(SVD))=0;           % ugly solution
SVD_col = squareform3d(SVD); % edge2column
SVD_col(SVD_col==0)=nan;     % ugly solution


%% add Yeo RNS labls
if not(isfield(cleanMat,'RNS'))
    regionsYeo = load(path_yeo);
    [LIA,LOCB] = ismember(cleanMat.regionDescriptions, regionsYeo.regionDescription);
    cleanMat.RNS = zeros(size(cleanMat.regionDescriptions));
    cleanMat.RNS(not(LIA)) = max(regionsYeo.RSN)+1;
    cleanMat.RNS(LIA) = regionsYeo.RSN(LOCB(LIA));
    
    % only when exist
    cleanMat.rsnDescription = [regionsYeo.rsnDescription; {'SubCortical'}; ...
        {'None'}];
end

% also subset
if not(isfield(cleanMat.NKE,'RNS'))
    [LIA,LOCB] = ismember(cleanMat.regionDescriptions, cleanMat.NKE.regionDescriptions);
    cleanMat.NKE.RNS = cleanMat.RNS(LIA);
end

%% Yeo-map of insomnia regions
cmRNS = ptc12(8, 'check');
cmRNS = cmRNS(1:end-1,:);  % ptc adds grey as n+1, will confuse plotBrain
    
% drop subcortical
regions=cleanMat.NKE.regionDescriptions;
values=cleanMat.NKE.RNS;

if strcmp(atlas, 'aparc')
    plotAtlas = 'aparc_aseg';
else
    plotAtlas = atlas;
end

% unavailable for Schaefer
% plotBrain(regions, values, cmRNS, ...
%     'atlas', plotAtlas, ...
%     'limits', [1,8], ...       %fixed to have similair map between atlasses
%     'Viewer', true, ...
%     'savePath', [result_path, '/yeo_insomnia_regions_colorblind']);



% (if any) set RNS nan values to 9 (corresponds to None Yeo label
cleanMat.NKE.RNS(isnan(cleanMat.NKE.RNS)) = 9;

% lausanne120 region index
region_number = 1:numel(cleanMat.NKE.RNS);

%  cluster (sort) RNS regions
[~,sort_index] = sort(cleanMat.NKE.RNS);

regionsTable = table(cleanMat.rsnDescription(cleanMat.NKE.RNS(sort_index)),...
    region_number', cleanMat.NKE.regionDescriptions(sort_index));

writetable(regionsTable,[result_path,'/regions_table.csv'],'Delimiter',',')  


%%% Edge estimates -----

% all labels are shuffled, p value calculated using permutations
% xMat with right side of the regression model
statsFA.regionsDescription = cleanMat.NKE.regionDescriptions;
statsFA.weight = {'FA'};

statsMD.regionsDescription = cleanMat.NKE.regionDescriptions;
statsMD.weight = {'MD'};

statsSVD.regionsDescription = cleanMat.NKE.regionDescriptions;
statsSVD.weight = {'SVD'};

%% GS - ID
% permutation is at 10, because we use tvalue instead of permuted p-value
statsFA.ID = permuteGSID(FA_col, xID, 10);  
statsMD.ID = permuteGSID(MD_col, xID, 10);
statsSVD.ID = permuteGSID(SVD_col, xID, 10);

%% GS - ITQ edge estimates
statsFA.ITQ = permuteGSITQ(FA_col, xITQ, 10);
statsMD.ITQ = permuteGSITQ(MD_col, xITQ, 10);
statsSVD.ITQ = permuteGSITQ(SVD_col, xITQ, 10);

%% ITQ edge specificity
% ITQ labels are shuffeled to determine how likely original t-values are due to random group
tic
statsFA.ITQperm = permuteITQ(FA_col, xITQ, nPerm);
statsMD.ITQperm = permuteITQ(MD_col, xITQ, nPerm);
statsSVD.ITQperm = permuteITQ(SVD_col, xITQ, nPerm);
toc

%% copy Yeo info
statsFA.RNS = cleanMat.NKE.RNS;
statsFA.rsnDescription = cleanMat.rsnDescription;

statsMD.RNS = cleanMat.NKE.RNS;
statsMD.rsnDescription = cleanMat.rsnDescription;

statsSVD.RNS = cleanMat.NKE.RNS;
statsSVD.rsnDescription = cleanMat.rsnDescription;

% save stats objects for reproducibility
save([result_path, '/', prefix, 'permutation_stats.mat'], ...
    'statsFA', 'statsMD', 'statsSVD', '-v7.3');

% to rerun a previous permutation
%load([result_path, '/', prefix, 'permutation_stats.mat'], 'statsFA', 'statsMD'); 



%% Build figures for the different modalities
modalities = {statsFA, statsMD, statsSVD};

for k = 1:length(modalities)
    statsObj = modalities{k};
    
    %% sign edge count and effect range
    t_values = cat(1, statsObj.ID.tVec(1,:), ...
                   squeeze(statsObj.ITQ.tVec(1,:,:))');
    sign_mask = abs(t_values) >= 2;
    beta_values = cat(1, statsObj.ID.estimateVec, statsObj.ITQ.estimateVec);

    T = sign_edge_count_range(beta_values, sign_mask, titles);
    writetable(T,[result_path, '/',statsObj.weight{:},'_', prefix, 'edge_table_'], 'FileType','spreadsheet');


    %% build circle plots with Yeo labeling
    tVals = cat(1, statsObj.ID.tVec(1,:), ...
                   squeeze(statsObj.ITQ.tVec(1,:,:))');
    bVals = cat(1, statsObj.ID.estimateVec, statsObj.ITQ.estimateVec);

     % find sign edges for all groups
     m = abs(tVals) >= 2;

     % plot sign betas
     sign_values = bVals .* m;
     [yeo_panel, ~] = yeo_schemaball_panel(sign_values, titles, cleanMat)

     % save as svg?
     print(yeo_panel, [result_path, '/',statsObj.weight{:},'_', prefix, 'structural_connectivity_', statsObj.weight{:}], ...
         '-dpng', '-r0');
     

    %% Functional profile barplots
    myYeo = unique(statsObj.RNS);
    tVals = cat(1, statsObj.ID.tVec(1,:), ...
                   squeeze(statsObj.ITQ.tVec(1,:,:))');
    bVals = cat(1, statsObj.ID.estimateVec, statsObj.ITQ.estimateVec);

    % find sign edges for all groups
    m = abs(tVals) >= 2;
    sign_values = bVals .* m;

    sign_matrix = squareform3d(sign_values);
    sign_regions = squeeze(any(sign_matrix, 2));

    % Yeo percentages
    for i = 1:length(myYeo)
        yeo_sign = sum(sign_regions(statsObj.RNS==myYeo(i),:), 1);
        yeo_total = sum(statsObj.RNS==myYeo(i));
        yeo_fa(i,:) = round(100 * yeo_sign / yeo_total);
    end

    % profile plots
    figure('color','w', 'Units','points' ,'Position', [5 5 2*539 360])
    yeo_labels = cleanMat.rsnDescription(myYeo);
    bar_colors = cmRNS(myYeo,:);
    for ii=1:6
            yeo_percentages = yeo_fa(:,ii);
            subplot(1,6,ii)

            b = bar(yeo_percentages,'facecolor', 'flat');
            set(gca, 'XTickLabel',yeo_labels, 'XTick',1:numel(yeo_labels));
            ylim([0,100])
            xtickangle(90)

            b.CData = bar_colors;
            split_title = split(titles(ii),',');
            title(split_title)
            set(gca,'box','off')
    end
    fig=gca;
    set(findall(fig, 'Type', 'Text'), 'FontSize', 12);
    print([result_path, '/',statsObj.weight{:},'_', prefix, 'func_profile'], ...
        '-dpng', '-r0');
end


% supplementary figures
for k = 1:length(modalities)
    statsObj = modalities{k};
    %% Number of deviating edges per functional network
    myYeo = unique(statsObj.RNS);
    tVals = cat(1, statsObj.ID.tVec(1,:), ...
                   squeeze(statsObj.ITQ.tVec(1,:,:))');

    % find sign edges for all groups
    sign_edges = abs(tVals) >= 2;
    sign_matrix = squareform3d(sign_edges);
    sign_edges_regions = squeeze(sum(sign_matrix, 2));

    % Yeo percentages
    for i = 1:length(myYeo)
        yeo_sign(i,:) = sum(sign_edges_regions(statsObj.RNS==myYeo(i),:), 1);
    end

    % abs edges profile plots
    figure('color','w', 'Units','points' ,'Position', [5 5 2*539 380])
    yeo_labels = statsObj.rsnDescription(myYeo);
    bar_colors = cmRNS(myYeo,:);
    yMax = max(yeo_sign,[],1:3);
    for ii=1:6
            yeo_edges = yeo_sign(:,ii);
            subplot(1,6,ii)

            b = bar(yeo_edges,'facecolor', 'flat');
            set(gca, 'XTickLabel',yeo_labels, 'XTick',1:numel(yeo_labels));
            ylim([0,yMax])
            xtickangle(90)

            b.CData = bar_colors;
            split_title = split(titles(ii),',');
            title(split_title)
            if ii==1
                ylabel('# of deviating connections')
            end
            set(gca,'box','off')
    end
    fig=gca;
    set(findall(fig, 'Type', 'Text'), 'FontSize', 10);
    print([result_path, '/',statsObj.weight{:},'_', prefix, 'edge_count_profile_suppl'], ...
        '-dpng', '-r0');
end


%% profile-based statistics
for k = 1:length(modalities)
    statsObj = modalities{k};
   


    myYeo = unique(statsObj.RNS);
    yeo_labels = statsObj.rsnDescription(myYeo);

    % create empty table
    profile_table = array2table(zeros(5,length(myYeo)+1), ...
        'VariableNames', [yeo_labels', 'profile'],...
        'RowNames', titles(2:end));


    % build Yeo distributions based ITQ permutations (only swaped itq labels)
    % loop over subtypes
    for i = 1:5
        perm_tVals = squeeze(statsObj.ITQperm.tVec(:,:,i));  % permutation vector
        m = abs(perm_tVals)>2;                              % edges with t value > 2
        sign_mat = squareform3d(m);                         % back to matrix
        sign_regions = squeeze(sum(sign_mat, 2)>0);         % simplify to if region contains sign edges

        % summarise to Yeo networks
        for ii = 1:length(myYeo)      
            % number of significant regions within Yeo network
            n_sign = sum(sign_regions(statsObj.RNS==myYeo(ii),:),1);
            yeo_hits(:,ii) = n_sign;        
        end

        % calculate p-values
        % per Yeo network (number of sign regions >= our result)
        yeo_prevalence = yeo_hits(2:end,:) >= yeo_hits(1,:);
        yeo_pval = sum(yeo_prevalence, 1) / nPerm;

        % odd of finding this profile
        profile_prevalence = all(yeo_prevalence,2);
        profile_pval = sum(profile_prevalence, 1) / nPerm;   
        profile_table{i,:} = [yeo_pval, profile_pval];

        clearvars yeo_hits
    end

    writetable(profile_table, ...
        [result_path, '/',statsObj.weight{:},'_', prefix, '_profile_statistics_'], ...
        'FileType','spreadsheet', 'WriteRowNames',true);
end
