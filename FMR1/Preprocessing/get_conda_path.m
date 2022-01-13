function out = get_conda_path()
    % Returns the path to conda.sh, path should be stored in conda_path.txt in
    % your matlab path

    out = fileread('conda_path.txt');
    out = strtrim(out);
