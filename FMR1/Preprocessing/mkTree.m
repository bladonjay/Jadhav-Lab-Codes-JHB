function out = mkTree(outDir)
    % mkTree(directory)
    % recursively make dirs so that the target folder is created

    if ~exist(fileparts(outDir),'dir')
        mkTree(fileparts(outDir));
    end
    out = 0;
    if ~exist(outDir,'dir')
        out = mkdir(outDir);
    else
        disp('Directory already exists')
    end
