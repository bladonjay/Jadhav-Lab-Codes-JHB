function out = ml_sort_on_segs(tetResDir,varargin)
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
    %   firings_out     : output filename for firings (default = 'tet/results/dir/firings_raw.mda')
    %   param_file      : parameter json file to use (default = 'tet/results/dir/params.json')
    %   time_file       : timestamps mda file for the day  (default = 'day/results/dir/animID_day_date.timestamps.mda')
    %   metrics_out     : output filename for cluster & isolation metrics (defualt = 'tet/results/dir/metrics_raw.json')
    %   tagged_metrics_out: output filename for curated cluster metrics (default = 'tet/res/dir/metrics_tagged.json')
    %   delete_temporary: flag whether to delete temporary metrics (uncombined cluster and isolation metrics) and cluster distance matrices from drift tracking (default = 1)
    %   firing_rate_thresh  : for curation, firing rate must be above this (default = 0.05)
    %   isolation_thresh    : for curation, isolation value must be above this (default = 0.95)
    %   noise_overlap_thresh: for curation, noise overlap must be below this (default = 0.03)
    %   peak_snr_thresh     : for curation, peak snr must be above this (default = 1.5)
    %   curation defaults are built into the ms4alg.create_label_map processor
    %   TODO: Cut into separate functions. Add overwrite flags

    if tetResDir(end)==filesep
        tetResDir = tetResDir(1:end-1);
    end

    geom = []; % optional csv defining electrode geometry (not needed for tetrodes)
    adjacency_radius = -1; % use all channels are one neighborhood (for tetrodes)
    detect_sign=-1; % sign of spikes to detect, Trodes 1.7+ automatically inverts sign on extraction so sign is +1
    detect_threshold = 3; % detection threshold for spike in st. dev from mean
    samplerate = 30000;

    % Curation parameters (from FrankLab)
    firing_rate_thresh = 0.01;
    isolation_thresh = 0.95;
    noise_overlap_thresh = 0.03;
    peak_snr_thresh = 1.5;

    % Automatic Filenames
    firings_out = [tetResDir filesep 'firings_raw.mda'];
    timeseries = [tetResDir filesep 'pre.mda'];
    if ~exist(timeseries,'file')
        timeseries = [tetResDir filesep 'pre.mda.prv'];
    end
    param_file = [tetResDir filesep 'params.json'];
    time_file = dir([fileparts(tetResDir) filesep '*timestamps*']);
    time_file = [time_file.folder filesep time_file.name];
    if strcmpi(time_file(end-2:end),'prv')
        time_prv =  jsondecode(fileread(time_file));
        time_file = time_prv.original_path;
    end
    metrics_out = [tetResDir filesep 'metrics_raw.json'];
    tagged_metrics_out = [tetResDir filesep 'metrics_tagged.json'];
    delete_temporary = 1;

    % Epoch detection parameters
    epoch_offsets = [];
    min_epoch_gap = 1; % seconds

    overwrite = 0;

    assignVars(varargin)

    if isempty(geom)
        geomStr = '';
    else
        geomStr = geom;
    end
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
    if exist(firings_out,'file') && ~overwrite
        fprintf('\n------\nFirings output already exists. Skipping...\n------\n')
    else
        % Determine epoch offsets (in samples) from timestamps.mda using gaps
        % larger than min_epoch_gap to separate epochs. Skip this is user provides
        % epoch offsets list (epoch_offsets)
        if isempty(epoch_offsets)
            timeDat = readmda(time_file);
            gaps = diff(timeDat);
            epoch_end = find(gaps>=min_epoch_gap*samplerate);
            epoch_offsets = [0 epoch_end];
            total_samples = numel(timeDat);
        else
            total_samples = numel(readmda(time_file));
        end

        fprintf('\n------\nProcessing %s\nDetected %i epochs. Extracting and Sorting\n------\n',tetResDir,numel(epoch_offsets));

        % Split into epoch segments and sort
        if numel(epoch_offsets)>1
            epoch_timeseries = cell(numel(epoch_offsets),1);
            epoch_firings = cell(numel(epoch_offsets),1);
            for k=1:numel(epoch_offsets)

                % Extract epoch timeseries
                t1 =  epoch_offsets(k);
                if k==numel(epoch_offsets)
                    t2 = total_samples-1;
                else
                    t2 = epoch_offsets(k+1)-1;
                end
                fprintf('\n------\nExtracting epoch %i/%i. Time: %0.2f s - %0.2f s\n------\n',k,numel(epoch_offsets),t1/samplerate,t2/samplerate);
                tmp_timeseries = sprintf('%s%spre-%02i.mda',tetResDir,filesep,k);
                extractInputs.timeseries = timeseries;
                extractOutputs.timeseries_out = tmp_timeseries;
                extractParams = struct('t1',t1,'t2',t2);
                ml_run_process('pyms.extract_timeseries',extractInputs,extractOutputs,extractParams);
                epoch_timeseries{k} = tmp_timeseries;

                % Sort epoch segment
                fprintf('\n------\nSorting epoch %i/%i. Paramters:\n\tadjacency_radius:%g\n\tdetect_sign:%i\n\tdetect_threshold:%g\n\tsamplerate:%g\n\tgeom: %s\n------\n',k,numel(epoch_offsets),adjacency_radius,detect_sign,detect_threshold,samplerate,geomStr);
                tmp_firings = sprintf('%s%sfirings-%02i.mda',tetResDir,filesep,k);
                sortInputs.timeseries = tmp_timeseries;
                if ~isempty(geom)
                    sortInputs.geom = geom;
                end
                sortOutputs.firings_out = tmp_firings;
                ml_run_process('ms4alg.sort',sortInputs,sortOutputs,sortParams);
                epoch_firings{k} = tmp_firings;
            end

            % anneal segments
            % TODO: Check epoch firings before running, errors if no spikes
            fprintf('\n------\nAnnealing Segments\n------\n')
            annealInputs.timeseries_list = epoch_timeseries;
            annealInputs.firings_list = epoch_firings;
            annealOutputs.firings_out = firings_out;
            offsetStr = sprintf('%i,',epoch_offsets);
            offsetStr = offsetStr(1:end-1);
            annealParams.time_offsets = offsetStr;
            % useless outputs, but required for process to run (errors if left empty)
            annealOutputs.dmatrix_out = [tetResDir filesep 'trash_dmatrix.mda'];
            annealOutputs.dmatrix_templates_out = [tetResDir filesep 'trash_dmatrix_templates.mda'];
            annealOutputs.k1_dmatrix_out = [tetResDir filesep 'trash_k1_dmatrix.mda'];
            annealOutputs.k2_dmatrix_out = [tetResDir filesep 'trash_k2_dmatrix.mda'];
            ml_run_process('pyms.anneal_segments',annealInputs,annealOutputs,annealParams);
        else 
            % if only 1 epoch
            fprintf('\n------\nSorting epoch. Paramters:\n\tadjacency_radius:%g\n\tdetect_sign:%i\n\tdetect_threshold:%g\n\tsamplerate:%g\n\tgeom: %s\n------\n',adjacency_radius,detect_sign,detect_threshold,samplerate,geomStr);
            sortInputs.timeseries = timeseries;
            sortOutputs.firings_out = firings_out;
            if ~isempty(geom)
                sortInputs.geom=geom;
            end
            ml_run_process('ms4alg.sort',sortInputs,sortOutputs,sortParams);
        end
        % firings output file has array NxL where the rows are
        % channel_detected_on,timestamp,cluster_labels and L is num data points
    end

    if  ~exist(metrics_out,'file') || overwrite
        % Compute cluster metrics
        fprintf('\n------\nComputing Cluster and Isolation Metrics\n------\n')
        metInputs.firings = firings_out;
        metInputs.timeseries = timeseries;
        metOutputs.cluster_metrics_out = [tetResDir filesep 'trash_metrics1.json'];
        isoOutputs.metrics_out = [tetResDir filesep 'trash_metrics2.json'];
        ml_run_process('ms3.cluster_metrics',metInputs,metOutputs,metParams);
        ml_run_process('ms3.isolation_metrics',metInputs,isoOutputs,isoParams);
        combineInput.metrics_list = {metOutputs.cluster_metrics_out;isoOutputs.metrics_out};
        combineOutput.metrics_out = metrics_out;
        ml_run_process('ms3.combine_cluster_metrics',combineInput,combineOutput)
    else
        fprintf('\n------\nCluster metrics already exist. Skipping creation...\n------\n')
    end
    % Add Curation Tags 
    % error in ms4alg.create_label_map: curation_spec.py.mp so skipping curation for now (9/13/18 RN)
    % Now using franklab's pyms.add_curation_tags
    if ~exist(tagged_metrics_out,'file') || overwrite
        fprintf('\n------\nAdding curation tags. Curation Parameters:\n\tfiring_rate_thresh: %f\n\tisolation_thresh\n\t %f\n\tnoise_overlap_thresh: %f\n\tpeak_snr_thresh: %f\n------\n',firing_rate_thresh,isolation_thresh,noise_overlap_thresh,peak_snr_thresh);
        pName = 'pyms.add_curation_tags';
        curInputs = struct('metrics',metrics_out);
        curOutputs = struct('metrics_tagged',tagged_metrics_out);
        ml_run_process(pName,curInputs,curOutputs,curParams);
    else
        fprintf('\n------\nAutomatically curated metrics already exist. Skipping curation...\n')
    end

    out = {firings_out;metrics_out;tagged_metrics_out};

    % delete temporary distance matrices
    if delete_temporary
        fprintf('\n------\nDeleting Extraneous Outputs: temporary metrics and anneal distance matrices\n------\n')
        delete([tetResDir filesep 'trash*'])
    end

    % delete epoch segment files
    if numel(epoch_offsets)>1
        fprintf('\n------\nDeleting temporary segment mda file\n----\n')
        for k=1:numel(epoch_offsets)
            delete(epoch_timeseries{k})
            delete(epoch_firings{k})
        end
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
