
function [statsObj] = permuteITQ(edgeMat, dataVec, nPerm)
%function [statsObj] = permuteITQ(edgeMat, groupVec, covarVec, nPerm)
% edgeMat should be 2d
n_groups = 5; % without control

% find edges
edge_index = find(sum(abs(edgeMat)>0, 1)>0);

% find controls
control_index = all(dataVec(:,1:n_groups)==0,2);
n_control = sum(control_index);
n_dsorder = sum(not(control_index));

% reorder data
var_control = dataVec(control_index,:);
var_ITQ = dataVec(not(control_index),:);
var_reordered = cat(1, var_control, var_ITQ);
edgeMat_reordered = cat(1, ...
    edgeMat(control_index,:),edgeMat(not(control_index),:));

% output structures
estimateVec = zeros(n_groups, length(edge_index));
tVec = zeros(nPerm+1, length(edge_index), n_groups);

% residual matrix and covar estimates
residMat = zeros(size(edgeMat));
bMat = zeros(size(dataVec,2)-n_groups, length(edge_index));

%% unpermuted effect
fprintf('Calculate effects \n')
for i = 1:length(edge_index)                       % across every edge
    y = edgeMat_reordered(:, edge_index(i));       % set edge values as y
    y(not(isnan(y))) = zscore(y(not(isnan(y))));   % zscore - ignore nan
    
    % full model
    mdl = regstats(y,var_reordered,'linear', 'tstat');
    estimateVec(:,i) = mdl.tstat.beta(2:1+n_groups);
    tVec(1,i,:) = mdl.tstat.t(2:1+n_groups);
    
    % reduced model    
    covmdl = regstats(y,var_reordered(:,1+n_groups:end),'linear', 'beta');
    b = covmdl.beta(2:end);
    bMat(:,edge_index(i)) = b;
    
    % residual without covar signal
    residMat(:, edge_index(i)) = y - (var_reordered(:,1+n_groups:end) * b);
end


%% split and reorder data
% enables permutation of only ITQ labels
edge_control = residMat(1:n_control,:);
edge_ITQ = residMat(1+n_control:end,:);
                                                                            
%% run permutations
fprintf('Permute ITQ labels \n')
% parallel - 6 workers
parfor (i = 1:nPerm, 6)
    fprintf(' permutation %d \n', i)
    
    % permute ITQ labels
    resid_edge = cat(1, edge_control, edge_ITQ(randperm(n_dsorder)',:));
    % Add permuted residual back to covar signal
    perm_edge = resid_edge + (var_reordered(:,1+n_groups:end)*bMat);    
    
    % edges loop
    tmpT = zeros(1, length(edge_index),n_groups); % tmp storage for loop output
    for ii = 1:length(edge_index)
        % set edge values as y
        y = perm_edge(:, edge_index(ii));
        
        % full model of permuted data
        mdl = regstats(y, var_reordered, 'linear', 'tstat');            
        tmpT(:,ii,:) = mdl.tstat.t(2:1+n_groups);            
        
    end
   % store edge level t-values for each permutation
   % (row 1 = original tvalues)
   tVec(1+i,:,:) = tmpT;
end

%% uncorrected p-values
% permutation based
abs_t = abs(tVec(1,:,:));
abs_t(isnan(abs_t))=0;  % set nan to zero preventing false positives
abs_tperm = abs(tVec(2:end,:,:));
p_perm = sum(abs_tperm >= abs_t, 1) / nPerm ;
p_perm = reshape(p_perm,length(edge_index),5)';
    
%% Romano-Wolf correction for multiple testing
%redundant
abs_tperm = abs(tVec(2:end,:,:));

% max t-value per permutation per ITQ
max_t_dist = max(abs_tperm,[],2);  

% compare original t-value to max T distribution
p_romano = sum(max_t_dist >= abs_t, 1) / nPerm ;
p_romano = reshape(p_romano,length(edge_index),5)';

%% return values to original structure
statsObj.estimateVec = zeros(n_groups, size(edgeMat,2));
statsObj.estimateVec(:,edge_index) = estimateVec;

statsObj.tVec = zeros(1+nPerm, size(edgeMat,2),n_groups);
statsObj.tVec(:,edge_index,:) = tVec;

statsObj.p_perm = zeros(n_groups, size(edgeMat,2));
statsObj.p_perm(:,edge_index) = p_perm;

statsObj.p_romano = zeros(n_groups, size(edgeMat,2));
statsObj.p_romano(:,edge_index) = p_romano;

statsObj.description = {'effect estimate of the linear regression'; ...
    'matrix with t-values for every edge in every permutation (row 1 is unpermuted)'; ...
    'p-value based on the permutations'; ...
    'p-value based on romano Wolf procedure (corrected for multiple comparison)'};
end