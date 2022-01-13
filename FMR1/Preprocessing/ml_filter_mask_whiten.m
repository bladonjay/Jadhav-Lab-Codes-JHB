function [out,maskingError] = ml_filter_mask_whiten(tetResDir,varargin)
    % outFile = ml_filter_mask_whiten(tetResDir,varargin) where tetResDir is the nt#.mountain
    % folder to be processed in the animal direct folder. should contain
    % raw.mda.prv and (optional) params.json
    % This function will use mountainLab bandpass filter the mda files, mask
    % out artifacts and whiten channels 
    % NAME-VALUE Pairs:
    %   freq_min    : min freq for bandpass filter (default = 300). can be overriden by params.json
    %   freq_max    : max freq for bandpass filter (default = 6000). can be overriden by params.json
    %   samplerate  : sampling rate of data in Hz  (default = 30000). can be overriden by params.json
    %   mask_artifacts: 1 to mask artifacts (default = 1)
    %   threshold   : number of st. deviations away from the mean to consider as artifacts (default=5). Default from Frank Lab
    %   chunk_size  : number of sample sets to zero if RSS is above threshold (default=2000). Default from Frank Lab

    freq_min = 600; % best for CA1
    freq_max = 6000;
    mask_artifacts = 1;
    samplerate = 30000;
    threshold = 5; % for artifact masking
    chunk_size = 2000; % for artifact masking

    filt_out = [tetResDir filesep 'filt.mda'];
    white_out = [tetResDir filesep 'pre.mda'];
    raw_timeseries = [tetResDir filesep 'raw.mda.prv'];
    overwrite = 0;

    assignVars(varargin)

    maskingError = 0;

    if tetResDir(end)==filesep
        tetResDir = tetResDir(1:end-1);
    end

    filtParams = struct('freq_min',freq_min,'freq_max',freq_max,'samplerate',samplerate);
    maskParams = struct('threshold',threshold,'chunk_size',chunk_size);

    if exist([tetResDir filesep 'params.json'],'file')
        paramTxt = fileread([tetResDir filesep 'params.json']);
        params = jsondecode(paramTxt);
        filtParams = setParams(filtParams,params);
        maskParams = setParams(maskParams,params);
    end

    if exist(filt_out,'file') && ~overwrite
        fprintf('\n------\nFiltered data already exists. Skipping filtering and masking...\n------\n')
    else
        % Bandpass filter
        fprintf('\n------\nFiltering, Maksing and Whitening for  %s\n\nBandpass Filtering. Parameters:\n\tfreq_min: %g\n\tfreq_max: %g\n\tsamplerate: %g\n------\n',tetResDir,freq_min,freq_max,samplerate)
        pName = 'ephys.bandpass_filter';
        inFile.timeseries = raw_timeseries;
        outFile.timeseries_out = filt_out;
        console_out = ml_run_process(pName,inFile,outFile,filtParams);

        % Mask Artifacts
        if mask_artifacts
            try
                fprintf('\n------\nMasking Artifacts. Parameters:\n\tthreshold: %g\n\tchunk_size: %i\n------\n',threshold,chunk_size);
                pName = 'ephys.mask_out_artifacts';
                inFile.timeseries = filt_out;
                outFile.timeseries_out = [tetResDir filesep 'filt_tmp.mda'];
                console_out = ml_run_process(pName,inFile,outFile,maskParams);
                delete(inFile.timeseries)
                movefile(outFile.timeseries_out,inFile.timeseries)
            catch ME
                fprintf('\n------\nMasking error: %s\n\n Skipping Masking...\n------\n',ME.message)
                maskingError = 1;
            end

        end
    end
    if exist(white_out,'file') && ~overwrite
        fprintf('\n------\nWhitened data already exists. Skipping whitening...\n------\n')
        out = white_out;
    else
        % Whiten
        fprintf('\n------\nWhitening Signal\n------\n');
        pName = 'ephys.whiten';
        inFile.timeseries = filt_out;
        outFile.timeseries_out = white_out;
        console_out = ml_run_process(pName,inFile,outFile);
        out = outFile.timeseries_out;
    end

function newParams = setParams(old,new)
    FNs = fieldnames(old);
    newParams = old;
    for k=1:numel(FNs)
        if isfield(new,FNs{k})
            newParams.(FNs{k}) = new.(FNs{k});
        end
        if isempty(newParams.(FNs{k}))
            newParams = rmfield(newParams,FNs{k});
        end
    end
    if isempty(fieldnames(newParams))
        newParams = [];
    end

