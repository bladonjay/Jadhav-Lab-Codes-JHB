function ml_process_animal(animID,rawDir,varargin)
    % ml_process_animal(animID,rawDir,varargin) will setup, preprocess and
    % spike sort data with mountainlab-js. rawDir is the raw data directory
    % containing the day directories. It is expected that mda files were
    % extracted from the rec files already.
    % Requires that day directories inside the rawDir are named day_date (e.g. 01_180819)
    % and that the mda directory inside the day dirs is labelled
    % animID_day_date.mda or animID_day_date_epoch.mda (the former will only
    % happen if you have 1 epoch), [e.g. RW9_02_1808224.mda]
    % This function will create a MountainSort directory in the direct folder.
    % Inside that directories will be .mountain directories for each day with
    % subdirectories for the sorting data for rode information will
    % be combined into a spikes mat-file in FilterFramework format for each
    % day. To avoid overwriting existing clustering, the spikes file will the
    % saved in the .mountain folder for each day.
    % NAME-VALUE Pairs:
    %   dataDir     : path to direct folder for animal. default = rawDir/../animID_direct
    %   sessionNums : array of days to process. default = [] proeach tetrode. Logs will be
    % created in each tetrode directory, and finally tetcesses all days in rawDir
    %   tet_list    : array of tetrodes to cluster. default = [] processes all tetrodes available
    %   mask_artifacts      : flag whether to mask artifacts before whitening.
    %                         default = 1, but the function errors rarely on some tetrodes (no idea
    %                         why yet) so you can choose not to do this process
    %                         *Actually, currently the script is setup to automatically 
    %                         skip artfact masking if it errors, this is mentioned at the end 
    %                         of the log file. so Check Your Logs!

    if rawDir(end)==filesep
        rawDir = rawDir(1:end-1);
    end
    % this is the data directory, it will be rootdir/animalID_direct just
    % beside your animal directory (rootdir/animalID)
    dataDir = [fileparts(rawDir) filesep animID '_direct'];
    sessionNums = [];
    tet_list = [];
    mask_artifacts = 1;

    assignVars(varargin);

    dayDirs = dir(rawDir);
    daysToProcess = zeros(numel(dayDirs),1);
    for k=numel(dayDirs):-1:1
        if strcmpi(dayDirs(k).name,'.') || strcmpi(dayDirs(k).name,'..')
            daysToProcess(k) = -1;
        else
            % parsing the day directory info: pat adheres to...
            % [two numbers] _ [many following numbers]
            % this parses day field into two numbers, and date into whatever
            % follows
            pat = '(?<day>[0-9]{2})_(?<date>\d*)';
            parsed = regexp(dayDirs(k).name,pat,'names');
            if isempty(parsed)
                disp(['Could not parse directory name: ' dayDirs(k).name])
                dn = input('What day is this (#)?  ');
            else
                dn = str2double(parsed.day);
            end
            daysToProcess(k) = dn; % record the day that was processed
        end
    end
    dayDirs(daysToProcess<0) = [];
    daysToProcess(daysToProcess<0) = [];

    % TODO: fix this, not restricting days properly
    % this is supposed to determine whether the days you entered match up
    % to the days that are parsed from the filefolder names
    if ~isempty(sessionNums)
        missing = setdiff(sessionNums,daysToProcess);
        [rmvDays,rmv] = setdiff(daysToProcess,sessionNums);
        dayDirs(rmv) = [];
        daysToProcess(rmv) = [];
        if ~isempty(missing)
            disp('Could not find data for days:')
            disp(missing)
        end
    end
    
    % this is the full contents of your animal directory
    dayDirs = strcat({dayDirs.folder},filesep,{dayDirs.name},filesep);

    disp('Processing raw data with mountain lab. Processing:')
    disp(dayDirs')
    if ~isempty(tet_list)
        disp('restricting to tetrodes:')
        disp(tet_list')
    end

    % Run mda_util, returns list of tetrode results directories
    %resDirs = mda_util(dayDirs,'tet_list',tet_list,'dataDir',dataDir);
    resDirs = mda_util(dayDirs,'tet_list',tet_list);
    dayResDirs = cell(numel(dayDirs),1);
    dayIdx = 1;
    maskErrors = zeros(numel(resDirs),1); 

    % For each day and tet sort spikes
    for k=1:numel(resDirs)
        rD = resDirs{k};
        diary([rD filesep 'ml_sorting.log'])
        fprintf('\n\n------\nBeginning analysis of %s\nDate: %s\n\nBandpass Filtering, Masking out artifacts and Whitening...\n------\n',rD,datestr(datetime('Now'))); 
        
        % filter mask and whiten
        % TODO: Add check to make sure there is data in the mda files, maybe up in mda_util
        [out,maskErrors(k)] = ml_filter_mask_whiten(rD,'mask_artifacts',mask_artifacts);
        % returns path to pre.mda.prv file

        fprintf('\n\n------\nPreprocessing of data done. Written to %s\n------\n',out)

        fprintf('\n------\nBeginning Sorting and curation...\n------\n')

        % Sort and curate
        out2 = ml_sort_on_segs(rD);
        % returns paths to firings_raw.mda, metrics_raw.json and firings_curated.mda
        %fprintf('\n\nSorting done. outputs saved at:\n    %s\n    %s\n    %s\n',out2{1},out2{2},out2{3})
        fprintf('\n\n------\nSorting done. outputs saved at:\n    %s\n    %s\n------\n',out2{1},out2{2})

        % Delete intermediate files
        %if ~keep_intermediates
        %    tmpDir = '/tmp/mountainlab-tmp/';
        %    disp('Removing intermediate processing mda files...')
        %    delete([tmpDir '*.mda'])
        %end
        %
        % Create matclust params file: Not Needed, only for trying to view spikes in matclust 5/8/19
        %fprintf('\n\n------\nCreating Matclust Params and Waves files\n------\n')
        %out3 = generateMatclustFromMountainSort(rD);
        %fprintf('\n\n------\nFile creation done. outputs saved at:\n    %s\n    %s\n------\n',out3{1},out3{2})

        if maskErrors(k)
            fprintf('\n######\nMasking error for this day. Masking Artifacts skipped. Spikes may be noisy\n######\n')
        end

        diary off

        % check if container results folders is already in dayResDirs and add if not
        if rD(end)==filesep
            rD = rD(1:end-1);
        end
        dD = fileparts(rD);
        if ~any(strcmpi(dayResDirs,dD))
            dayResDirs{dayIdx} = dD;
            dayIdx = dayIdx+1;
        end
    end
    fprintf('Completed automated clustering!\n')
    
    %% Converting firings to spikes FF shoould be done after manual cluster verification with qt-mountainview
    %for k=1:numel(dayResDirs)
    %    dD = dayResDirs{k};
    %    fprintf('Creating spikes file for %s...\n',dD);
    %    [remainder,dirName] = fileparts(dD);
    %    if isempty(dirName)
    %        [~,dirName] = fileparts(remainder);
    %    end
    %    pat = '\w*_(?<day>[0-9]+)_\w*';
    %    parsed = regexp(dirName,pat,'names');
    %    dayNum = str2double(parsed.day);
    %    fprintf('Identified Day Number as %02i\n',dayNum);
    %    spikesFile = convert_ml_to_FF(animID,dD,dayNum);
    %    fprintf('Done! Created %s\n\n',spikesFile)
    %end








        

