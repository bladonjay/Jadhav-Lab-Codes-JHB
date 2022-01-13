function out = ml_sort_full(tetResDir,varargin)
    % out = ml_sort_on_segs(tetResDir,varargin) will spike sort the data in the
    % given folder where tetResDir is the .nt#.mountain folder in the direct
    % directory and contains pre.mda.prv and (optional) params.json
    % This function will sort spikes, compute cluster metrics and create
    % curation tags. All params can be oveerriden via Name-Value pairs in the
    % input or via the params.json file (which has priority)
    % NAME-VALUE Pairs:
    %   geom    : path to csv file with electrode geometry (default = [] for tetrodes)
    %   adjacency_radius    : see ml-spec ephys.ms4alg -p (default = -1 for tetrodes)
    %   detect_sign     : sign of spikes to detect (default = 1)
    %   samplerate      : sampling rate of the data in Hz (default = 30000)
    %   firing_rate_thresh  : for curation, firing rate must be above this (default = 0.05)
    %   isolation_thresh    : for curation, isolation value must be above this (default = 0.95)
    %   noise_overlap_thresh: for curation, noise overlap must be below this (default = 0.03)
    %   peak_snr_thresh     : for curation, peak snr must be above this (default = 1.5)
    %   curation defaults are built into the ms4alg.create_label_map processor


    geom = []; % optional csv defining electrode geometry (not needed for tetrodes)
    adjacency_radius = -1; % use all channels are one neighborhood (for tetrodes)
    detect_sign=1; % sign of spikes to detect, Trodes 1.7+ automatically inverts sign on extraction so sign is +1
    detect_threshold = 3; % detection threshold for spike in st. dev from mean
    samplerate = 30000;
    % Curation parameters
    firing_rate_thresh = [];
    isolation_thresh = [];
    noise_overlap_thresh = [];
    peak_snr_thresh = [];
    firings_out = [tetResDir filesep 'firings_raw.mda'];
    metrics_out = [tetResDir filesep 'metrics_raw.json'];
    timeseries = [tetResDir filesep 'pre.mda'];

    assignVars(varargin)

    % Set parameters for sorting, metrics, and curation. Replace any parameters with contents of params.json
    sortParams = struct('adjacency_radius',adjacency_radius,'detect_sign',detect_sign,'detect_threshold',detect_threshold);
    metParams = struct('samplerate',samplerate);
    isoParams = struct('compute_bursting_parents','true');
    curParams = struct('firing_rate_thresh',firing_rate_thresh,'isolation_thresh',isolation_thresh,'noise_overlap_thresh',noise_overlap_thresh,'peak_snr_thresh',peak_snr_thresh);
    if exist(param_file,'file')
        paramTxt = fileread(param_file);
        params = jsondecode(paramTxt);
        if isfield(params,'samplerate')
            samplerate = params.samplerate;
        end
        sortParams = setParams(sortParams,params);
        metParams = setParams(metParams,params);
        isoParams = setParams(isoParams,params);
        curParams = setParams(curParams,params);
    end

    % Sort entire file at once

    % Sort 
    pName = 'ms4alg.sort';
    sortInputs.timeseries = timeseries;
    sortOutputs.firings_out = firings_out;
    if ~isempty(geom)
        sortInputs.geom = geom;
    end
    ml_run_process(pName,sortInputs,sortOutputs,sortParams);
    % output file have array NxL where the rows are
    % channels_used,timestamp,cluster_labels and L is num data points

    % Compute cluster metrics
    metInputs.firings = firings_out;
    metInputs.timeseries = timeseries;
    metOutputs.cluster_metrics_out = [tetResDir filesep 'trash_metrics1.json'];
    isoOutputs.metrics_out = [tetResDir filesep 'trash_metrics2.json'];
    ml_run_process('ms3.cluster_metrics',metInputs,metOutputs,metParams);
    ml_run_process('ms3.isolation_metrics',metInputs,isoOutputs,isoParams);
    combineInput.metrics_list = {metOutputs.cluster_metrics_out;isoOutputs.metrics_out};
    combineOutput.metrics_out = metrics_out;
    ml_run_process('ms3.combine_cluster_metrics',combineInput,combineOutput)

    % Add Curation Tags 
    % error in ms4alg.create_label_map: curation_spec.py.mp so skipping curation for now (9/13/18 RN)
    % Now using franklab's pyms.add_curation_tags
    pName = 'pyms.add_curation_tags';
    curInputs = struct('metrics',metrics_out);
    curOutputs = struct('metrics_tagged',metrics_out);
    ml_run_process(pName,curInputs,curOutputs,curParams);

    out = {sortOutputs.firings_out;metOutputs.metrics_out};

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
