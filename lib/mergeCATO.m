

function mergedMatrix = mergeCATO(catoDir, matPattern)
% function merges all .mat cato matrices in the specified dir into a single
% object based on the specified pattern. Files names are parse at '_' and
% the first element is used as subject name. 
%
% e.g. mergeCATO(mydir, *connectivity_csd_dti_lausanne120.mat) gathers all 
% csd_dti_lausanne120.mat matrices from mydir and combines them into a
% single object. _

% arguments
switch nargin
    case 0
        error('Error: Not enough arguments.');
    case 1
        if exist('catoDir', 'var')
            error('Error: Specify .mat pattern.')
        else
            catoDir = pwd ;
        end         
    case 2
        if not(isfolder(catoDir))
            error('Error: Cannot fir folder')
        end
    otherwise
        error('Error: Too many arguments.')
end

% list available .mat files
filePattern = fullfile(catoDir, matPattern);
theFiles = dir(filePattern);

% pre allocate output
mergedMatrix = struct();

    % loop
    for k = 1 : length(theFiles)
        % set file name and path
        baseFileName = theFiles(k).name;
        fullFileName = fullfile(theFiles(k).folder, baseFileName);
        fprintf(1, 'Reading %s\n', baseFileName);
        % load .mat file
        catoArray = load(fullFileName);

        % first itteration adds region descriptions etc
        if k == 1
            mergedMatrix.regionDescriptions = catoArray.regionDescriptions ;
            mergedMatrix.ROIs = catoArray.ROIs ;
            mergedMatrix.weightDescriptions = catoArray.weightDescriptions ;
            mergedMatrix.subjects = cell(length(theFiles), 1);
            mergedMatrix.connectivity = zeros(size(catoArray.connectivity, 1), ...
                size(catoArray.connectivity, 2), size(catoArray.connectivity, 3), ...
                length(theFiles));
        end

        % add connectivity matrices to correct field
        A = strsplit(baseFileName, '_');
        mergedMatrix.subjects(k,1) = A(1);
        mergedMatrix.connectivity(:,:,:,k) = catoArray.connectivity;
    end
end




