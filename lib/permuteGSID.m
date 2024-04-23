
function [statsObj] = permuteGSID(edgeMat, dataVec, nPerm)
% edgeMat should be 2d
n_groups = 1; % without control

% find edges
edge_index = find(sum(abs(edgeMat)>0, 1)>0);

% output structures
estimateVec = zeros(1, length(edge_index));
tVec = zeros(nPerm+1, length(edge_index));
p_lm = zeros(1, length(edge_index));

% residual matrix and covar estimates
residMat = zeros(size(edgeMat));
bMat = zeros(size(dataVec,2)-n_groups, length(edge_index));

% unpermuted effect
fprintf('Calculate effects \n')
for i = 1:length(edge_index)                       % across every edge
    y = edgeMat(:, edge_index(i));                 % set edge values as y
    y(not(isnan(y))) = zscore(y(not(isnan(y))));   % zscore - ignore nan
    
    % full model
    mdl = regstats(y,dataVec,'linear', 'tstat');
    estimateVec(:,i) = mdl.tstat.beta(2:1+n_groups);
    tVec(1,i) = mdl.tstat.t(2:1+n_groups);
    p_lm(:,i) = mdl.tstat.pval(2:1+n_groups);
    
    % reduced model    
    covmdl = regstats(y,dataVec(:,1+n_groups:end),'linear', 'beta');
    b = covmdl.beta(2:end);
    bMat(:,edge_index(i)) = b;
    
    % residual without covar signal
    residMat(:, edge_index(i)) = y - (dataVec(:,1+n_groups:end) * b);   
end

                                                                               
%% run permutations
fprintf('Permute all labels \n')
n_all = size(dataVec,1);

% parallel - 6 workers
parfor (i = 1:nPerm, 6)
    fprintf(' permutation %d \n', i)
    
    % permute ITQ labels
    perm_resid = residMat(randperm(n_all)',:);
    % Add permuted residual back to covar signal
    perm_edges = perm_resid + (dataVec(:,1+n_groups:end)*bMat);
    
    % edges loop
    tmpT = zeros(1, length(edge_index)); % tmp storage for loop output
    
    for ii = 1:length(edge_index)
        y = perm_edges(:, edge_index(ii));      % set edge values as y  
        mdl = regstats(y, dataVec, 'linear', 'tstat');            
        tmpT(:,ii) = mdl.tstat.t(2:1+n_groups); 
    end
   % store edge level t-values for each permutation
   % (row 1 = original tvalues)
   tVec(1+i,:) = tmpT;
end

%% uncorrected p-values
% permutation based
abs_t = abs(tVec(1,:));
abs_t(isnan(abs_t))=0; % set nan to zero preventing false positives
abs_tperm = abs(tVec(2:end,:));
p_perm = sum(abs_tperm >= abs_t, 1) / nPerm ;
    
%% Romano-Wolf correction for multiple testing
%redundant
abs_tperm = abs(tVec(2:end,:,:));

% max t-value per permutation per ITQ
max_t_dist = max(abs_tperm,[],2);  

% compare original t-value to max T distribution
p_romano = sum(max_t_dist >= abs_t, 1) / nPerm ;

statsObj.estimateVec = zeros(n_groups, size(edgeMat,2));
statsObj.estimateVec(:,edge_index) = estimateVec;

statsObj.tVec = zeros(1+nPerm, size(edgeMat,2),n_groups);
statsObj.tVec(:,edge_index,:) = tVec;

statsObj.p_lm = zeros(n_groups, size(edgeMat,2));
statsObj.p_lm(:,edge_index) = p_lm;

statsObj.p_perm = zeros(n_groups, size(edgeMat,2));
statsObj.p_perm(:,edge_index) = p_perm;

statsObj.p_romano = zeros(n_groups, size(edgeMat,2));
statsObj.p_romano(:,edge_index) = p_romano;

statsObj.description = {'effect estimate of the linear regression'; ...
    'matrix with t-values for every edge in every permutation (row 1 is unpermuted)'; ...
    'p-values of the linear regression model'; ...
    'p-value based on the permutations'; ...
    'p-value based on romano Wolf procedure (corrected for multiple comparison)'};
end