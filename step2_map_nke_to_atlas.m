%% env
% set variables
projdir = '/localpath';
atlas = 'lausanne120'; % options: lausanne120 aparc schaefer100-yeo7
idClusterPath = [projdir, '/data/insomnia_clusters_k25-50.csv'];
figurepath = [projdir,'/figures'];

% missing parcellation_ATLAS files can be generated using make_parcellation_ATLAS.sh
nke_data_path = '/localpath/nke-viewer/data';
datapath = '/localpath/parcellation/parcellation_ATLAS.nii.gz';
examplepath = '/localpath/parcellation/example_connectivity_dti_ATLAS.mat';

% add paths
addpath(genpath(projdir));
addpath(genpath('/localpath/Simple-Brain-Plot'));


% set paths
datapath = strrep(datapath, 'ATLAS', atlas);
examplepath = strrep(examplepath, 'ATLAS', atlas);



%% load 
% atlas ROI numbers and labels
atlasFeatures = load(examplepath, 'ROIs', 'regionDescriptions');

% load atlas in MNI152 space
ref = load_nifti(datapath);
ref.vol = double(ref.vol); % all values must be same type for accumarray.

% insomnia clusters per cluster_grouping
allCOI = readtable(idClusterPath);



%% loop over all cluster groupings

% empty object to store insomnia clusters of each cluster grouping
heatmapData = zeros(length(atlasFeatures.ROIs), length(allCOI.cluster_grouping));

% loop over every cluster_grouping
for k=1:length(allCOI.cluster_grouping)
    % set nfeatures as cluster_grouping
    nfeatures=allCOI.cluster_grouping{k}; 
    
    % set path to data files of specific cluster_grouping
    feature_path = fullfile(nke_data_path, sprintf('%s', nfeatures));
    
    % extract row from table with clusters per cluster_grouping
    current = allCOI(strcmp(allCOI.cluster_grouping, sprintf('%s', nfeatures)), :);
    COI = str2num(current.clusters{:});

    % empty vector for output
    ROIs_per_cluster = zeros(length(atlasFeatures.ROIs), length(COI));

    % loop over clusters of interest
    for ic = COI
        % load NKE cluster of interest in MNI152 space
        feature = load_nifti(fullfile(feature_path, sprintf('circuit_%s_dom%.2i.nii.gz', nfeatures, ic)));
        feature.vol = double(feature.vol);

        % create mask
        feature.mask = feature.vol > 0 ;

        % mask atlas refference (ROI labels)
        filtered = ref.vol .* feature.mask ;

        % store ROIs associated with cluster/domain
        ROIs_per_cluster(:, COI == ic) = ismember(atlasFeatures.ROIs, unique(filtered)) ;
    end



    %% Plot (not dynamic with atlas variable yet)
    load('regionDescriptions.mat');

    % sum the columns to label regions that are involved in multiple
    values = sum(ROIs_per_cluster, 2);

    % add values to heatmap data
    heatmapData(:, k) = values;
    
    % if lausanna120, match atlas regions using labels 
        %lausanna120 has less regions than lausanne120_aseg in plotBrain
    if strcmp(atlas, 'lausanne120')
        % cross with aparc_aseg to drop aseg regions
        ind = ~ismember(regionDescriptions.lausanne120_aseg, regionDescriptions.aparc_aseg);
        regions = regionDescriptions.lausanne120_aseg(ind) ;
    else
        regions = atlasFeatures.regionDescriptions;
    end
        
    % set color scheme
    cm = cbrewer('qual', 'Set3', 5);
    cm(1,:) = [1    1    1]; %white

    % test should be same length as rDescr with  a value  for the regions I
    % want coloured
    % plotBrain(regions, values, cm, ...
    % 'atlas', [atlas], ...
    % 'Viewer', false, ...
    % 'savePath', sprintf('/Users/tombresser/Documents/PhD/proj-subtypes/figures/%s/%s', atlas, nfeatures));
end

%% visualize heatplot off all k mappings (overlapping)

values = sum(heatmapData, 2);
cm = cbrewer('seq', 'YlOrRd', max(values));
cm(1,:) = [1    1    1]; %white
% 
% plotBrain(regions, values, cm, ...
%     'atlas', [atlas], ...
%     'Viewer', true, ...
%     'savePath', [figurepath,'/', atlas,'__k25-50_heatmap_all']);


%% Make heatmap figure showing regions >= 10 insomnia domains
insomniaROIs = struc([]);
insomniaROIs.regions = regions;
insomniaROIs.values = sum(heatmapData, 2);

regions = insomniaROIs.regions;

% select regions that link at least 10 times to any neurobiological domains
%   associated with insomnia subtype traits
n = 10;
values = insomniaROIs.values;
values(values < n) = 0;

plotBrain(regions, values, cm, ...
    'atlas', atlas, ...
    'Viewer', true, ...
    'savePath', [figurepath,'/', atlas,'__k25-50_heatmap']);
    
regionList = regions(values > 0);
    
fileID = fopen([projdir,'/data/',atlas,'_regions_threshold',num2str(n),'.txt'],'w');
fprintf(fileID,'%s\n', regionList{:});
fclose(fileID);


