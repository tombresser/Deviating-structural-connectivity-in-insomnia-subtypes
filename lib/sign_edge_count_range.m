
function [T] = sign_edge_count_range(beta_values, sign_mask, titles)
% initiate table
headers = {'Subtype' 'posEdge' 'posRange' 'negEdge' 'negRange'};
data = cell(length(titles),5);
T = cell2table(data);
T.Properties.VariableNames = headers;
T.Subtype = titles';   


% create panel
for j = 1:length(titles)
    sign_beta = beta_values(j,:) .* sign_mask(j,:);

    % table values
    posEdge = sign_beta(sign_beta > 0);
    T.posEdge(j) = {numel(posEdge)};
    [S,L] = bounds(posEdge);
    T.posRange(j) = {sprintf('%0.3f to %0.3f', S, L)};

    negEdge = sign_beta(sign_beta < 0);
    T.negEdge(j) = {numel(negEdge)};
    [S,L] = bounds(negEdge);
    T.negRange(j) = {sprintf('%0.3f to %0.3f', S, L)};

    clearvars S L sign_beta posEdge negEdge
end