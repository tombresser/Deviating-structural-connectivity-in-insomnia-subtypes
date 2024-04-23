function [abbr_names] = abbreviate_DK(region_names)
    % replaces regions names with abbreviations assuming Desikan-Killiany
    % naming. Also works for the extended lausanne parcellations.
    % 
    % Abbreviations used as found in doi.org/10.3390/brainsci11010107
    % https://www.researchgate.net/figure/Regions-in-the-Desikan-Killiany-DK-atlas_tbl2_348518431
    % 
    % Two additions used in the larger lausanne parcellations:
    % 'unknown' = 'Unkn'
    % 'corpuscallosum' = 'CC'
    
    atlas_pattrn = {'bankssts', 'caudalanteriorcingulate', ...
        'caudalmiddlefrontal', 'cuneus', 'entorhinal', 'frontalpole', ...
        'fusiform', 'inferiorparietal', 'inferiortemporal', 'insula', ...
        'isthmuscingulate', 'lateraloccipital', 'lateralorbitofrontal', ...
        'lingual', 'medialorbitofrontal', 'middletemporal', 'paracentral', ...
        'parahippocampal', 'parsopercularis', 'parsorbitalis',  ...
        'parstriangularis', 'pericalcarine', 'postcentral', ...
        'posteriorcingulate', 'precentral', 'precuneus', ...
        'rostralanteriorcingulate', 'rostralmiddlefrontal', 'superiorfrontal', ...
        'superiorparietal', 'superiortemporal', 'supramarginal', ...
        'temporalpole', 'transversetemporal', 'unknown', 'corpuscallosum'};
    
    abbrev = {'B', 'CACg', 'CMF', 'Cu', 'En', 'FPol', 'Fu', 'IP', 'IT', 'Ins', 'IstCg', ...
    'LO', 'LOrF', 'Lg', 'MOrF', 'MT', 'PaC', 'PaH', 'Op', 'Or', 'Tr', ...
    'PerCa', 'PoC', 'PoCg', 'PreC', 'PreCu', 'RoACg', 'RoMF', 'SF', ...
    'SP', 'ST', 'SM', 'TPol', 'TrT', 'Unkn', 'CC'};
    
    % replace core name with abbreviation
    abbr_names = replace(region_names, atlas_pattrn, abbrev);
    
    % simplify cortex + left/right with one letter
    abbr_names = strrep(abbr_names, 'ctx-rh-', 'r-');
    abbr_names = strrep(abbr_names, 'ctx-lh-', 'l-');
    abbr_names = strrep(abbr_names, 'rh-', 'r-');
    abbr_names = strrep(abbr_names, 'lh-', 'l-');
    abbr_names = strrep(abbr_names, 'rh_', 'r-');
    abbr_names = strrep(abbr_names, 'lh_', 'l-');
    
    % drop addional _'s (used when regions are subdivided e.g. _1 and _2)
    abbr_names = strrep(abbr_names, '_', '');    




