function [fig, cBar1, cBar2] = yeo_schemaball_panel(values, titles, roiStruc, varargin)
% assumes effect sizes within -1 to 1

% Default values
wholebrainFlag = false;

% Check if 'wholebrain' flag is provided
if ~isempty(varargin)
    % Loop through varargin to find switches
    for i = 1:numel(varargin)
        if strcmpi(varargin{i}, 'wholebrain')
            wholebrainFlag = true;
        end
    end
end 

% set vars
yeoLabels=roiStruc.rsnDescription;
if wholebrainFlag == true
    RNS=roiStruc.RNS;
    regions=roiStruc.regionDescriptions;
else
    RNS=roiStruc.NKE.RNS;
    regions=roiStruc.NKE.regionDescriptions;
end

% hard coded load
addpath('/local_path/subplot_tight')

% abbreviated region names
short_region_names = abbreviate_DK(regions);
short_region_names = pad(short_region_names, max(strlength(short_region_names)));

% try to cluster (sort) RNS regions
[~,sort_index] = sort(RNS);

%% color maps
% effect size
cm=cbrewer('div', 'RdBu', 1000); 
cm=flipud(cm);

% Colorscheme based on RSN ( Yeo )
cmRNS = ptc12(8);  % hardcoded to be similair across atlasses
cmYeo = cmRNS(RNS,:);

%% start figure

fig = figure('color', 'white', ...
    'Units','points' ,'Position', [5 5 2*539 2*360]);

for i = 1:length(titles)

    % matrix with edge values between regions
    W = squareform3d(values(i,:));
    subplot_tight(2,3,i,[0.001])
    
    schemaball(cellstr(string(1:size(W,1))), ...
       W(sort_index,sort_index), cm, -1, 1, cmYeo(sort_index,:))
    
    split_title = split(titles(i),',');
    title(split_title, 'Units', 'normalized', 'Position', [0.5, 0.9, 0]) 
end

% adjust fontsize
set(findall(fig, 'Type', 'Text'), 'FontSize', 10);


%% seperate color bars
% RNS
cmTMP = cmRNS(unique(RNS),:);
cBar1 = figure('color', 'white')
colormap(cmTMP)

% legend labels
rnsLabels = yeoLabels(unique(RNS))

colorbar('Ticks', linspace(0+0.1,1-0.1, length(cmTMP)), 'TickLabels', rnsLabels, ...
    'TickLength', 0, 'Position', [0.7 0.1 0.05 0.8], ...
    'FontSize', 12)
set(gca,'XColor', 'none','YColor','none')

% connectivity strength
cBar2 = figure('color', 'white')
colormap(cm)
colorbar('Ticks', linspace(0,1,11), ...
    'TickLabels', -1:0.2:1, ...
    'FontSize', 12)
set(gca,'XColor', 'none','YColor','none')
