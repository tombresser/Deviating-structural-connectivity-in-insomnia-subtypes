function M = neuroCombat_wrap(dataMat, modelMat, rPath, cachePath)
    % wrapper to run neuroCombat in R. Messges, warnings  and errors are
    % logged in projPath/log/. neuroCombat expects subjects as COLUMNs, but
    % the R script dubble checks this and pivots the data if needed. Output
    % is in the same dimension.

    % temp save model matrix
    writetable(modelMat, [cachePath '/tmp_combat_model.csv'])  
    
    % temp save data matrix
    writematrix(dataMat, [cachePath '/tmp_combat_data.csv']) 

    % run Rscript
    % expects tmp_combat_data.csv & tmp_combat_model
    system([rPath,' /some_path/GitDir/lib/ComBat_edges_matlab.r'])
    
    % read output
    M = readmatrix([cachePath '/tmp_combat_harmonized.csv']);
    
    % remove temp files
    delete([cachePath '/tmp_combat_model.csv'], ...
        [cachePath '/tmp_combat_data.csv'], ...
        [cachePath '/tmp_combat_harmonized.csv'])

