function out = convert_ml_to_FF_nocuration(animID,animDir,sessionNums,varargin)
%animDir is animal_direct folder with processed data.

%sessionNum is day number

samplerate = 30000;
min_epoch_break=1;
overwrite = 1;

assignVars(varargin)

%% Find correct day folder inside MountainSort results
dayfolders = dir([animDir filesep 'MountainSort']);
dayfolders = {dayfolders(3:end).name};

for s = 1:length(sessionNums)
    sessionNum = sessionNums(s);
    sessionstr = sprintf('%s_%02i',animID,sessionNum);
    resDir = [animDir, filesep, 'MountainSort', filesep, dayfolders{contains(dayfolders,sessionstr)}];
    
    
    %% Collect Data and create save file
    out = '';
    tetDirs = dir([resDir filesep '*.mountain']);
    timeFile = dir([resDir filesep '*timestamps*']);
    if isempty(timeFile)
        error('No timestamps file found in %s',resDir);
    end
    timeFile = [timeFile.folder filesep timeFile.name];
    if strcmp(timeFile(end-2:end),'prv')
        try
            tmp = jsondecode(fileread(timeFile));
            timeFile = tmp.original_path;
            timeDat = readmda(timeFile);
        catch
            %in case raw data is in a different location
            rawDir = '/media/jadhav/JDS_CAUSAL_DRIVE/PFC_SWRcausal/JS50/';
            daystr = erase(sessionstr,[animID,'_']);
            datestr = erase(erase(dayfolders{contains(dayfolders,sessionstr)}, [animID,'_',daystr,'_']),'.mountain');
            rawDir = [rawDir, animID,'/',daystr,'_',datestr,'/',animID,'_',datestr,'.mda/'];
            timeFile = dir([rawDir filesep '*timestamps*']);
            timeDat = readmda([rawDir filesep timeFile.name]);
        end
    end
    
    saveFile = sprintf('%s%s%sspikes%02i.mat',animDir,filesep,animID,sessionNum);
    if exist(saveFile,'file') && ~overwrite
        fprintf('Spikes file @ %s already exists. Skipping...\n',saveFile)
        return;
    end
    
    %% Determine epoch start times by looking for gaps longer than 1 sec
    gaps = diff(timeDat);
    epoch_gaps = find(gaps>=min_epoch_break*samplerate);
    epoch_starts = timeDat([1 epoch_gaps]);
    Nepochs = numel(epoch_starts);
    fprintf('Detected %i epochs. Start times:\n',Nepochs);
    disp(epoch_starts');
    
    %% get tet numbers
    pat = '\w*.nt(?<tet>[0-9]+).\w*';
    parsed = cellfun(@(x) regexp(x,pat,'names'),{tetDirs.name});
    tet_nums = str2double({parsed.tet});
    Ntets = max(tet_nums);
    fprintf('Detected %i tetrodes. Tetrodes:\n',Ntets);
    disp(tet_nums');
    drawnow;
    
    spikes = cell(1,sessionNum);
    spikes{sessionNum} = cell(1,Nepochs);
    spikes{sessionNum}(:) = {cell(1,Ntets)};
    
    %% Convert from firings data
    for k=1:numel(tetDirs)
        tD = [tetDirs(k).folder filesep tetDirs(k).name];
        tetNum = tet_nums(k);
        metFile = [tD filesep 'metrics_tagged.json'];
        
        %tetrodes with no cells will not have firings.curated.mda files
        if isfile([tD filesep 'firings_raw.mda'])
            fireFile = [tD filesep 'firings_raw.mda'];
            metDat = jsondecode(fileread(metFile));
            metDat = metDat.clusters;
            
            %CS- remove cells that are tagged as "rejected." Keep all other tags
            %and save in spikes file.
            tags = {metDat.tags};
            multiunit = cellfun(@(x) strcmp('mua',x),tags,'UniformOutput',false);
            acceptcells = ~cellfun(@any, multiunit);
            metDat = metDat(acceptcells);
            
            fireDat = readmda(fireFile); % Rows are # Channels, timestamp (starting @ 0), cluster #
            
            %if the first timestamp entry is 0, replace with 1
            
            if fireDat(2,1) == 0
                fireDat(2,1) = 1;
            end
            
            clusters = [metDat.label];
            %Nclust = max(clusters)+any(clusters==0);
            [clusters,ic] = sort(clusters,'ascend');
            metDat = metDat(ic);
            
            % TODO: check params.json for smaplerate and error if not matching
            
            fireTimes = timeDat(fireDat(2,:));
            for l=1:Nepochs
                spikes{sessionNum}{l}{tetNum} = cell(1,length(clusters));
                for m=1:length(clusters)
                    t1 = epoch_starts(l);
                    if l<Nepochs
                        t2 = epoch_starts(l+1)-1;
                    else
                        t2=timeDat(end);
                    end
                    
                    
                    idx = fireDat(3,:)==clusters(m) & fireTimes>=t1 & fireTimes<=t2;
                    mD = metDat(m).metrics;
                    dat = zeros(sum(idx),7);
                    fields = 'time x y dir not_used amplitude(highest variance channel) posindex n_detection_channels';
                    info = 'spike data from MountainSort 4 (MountainLab-JS)';
                    meanrate = mD.firing_rate;
                    isolation = mD.isolation;
                    num_events = mD.num_events;
                    noise_overlap =mD.noise_overlap;
                    peak_noise = mD.peak_noise;
                    peak_amplitude = mD.peak_amp;
                    peak_snr = mD.peak_snr;
                    %refractory_violation_1msec = mD.refractory_violation_1msec;
                    dat(:,1) = fireTimes(idx)'/samplerate;
                    dat(:,7) = fireDat(1,idx)';
                    %% add tag
                    
                    if isempty(dat)
                        spikes{sessionNum}{l}{tetNum}{m} = [];
                    else
                        
                        %                     tag = metDat(m).tags;
                        %                     if length(tag) > 1
                        %                         tag = tag{(~(cellfun('isempty',strfind(tag,'accepted'))))};
                        %                     end
                        
                        %populate fields
                        spikes{sessionNum}{l}{tetNum}{m} = struct('data',dat,'meanrate',meanrate,...
                            'info',info,'fields',fields,...
                            'timerange',[t1 t2]/samplerate,...
                            'isolation',isolation,...
                            'num_events',num_events,'noise_overlap',noise_overlap,...
                            'peak_amplitude',peak_amplitude,'peak_noise',peak_noise,...
                            'peak_snr',peak_snr);
                    end
                    
                end
            end
        end
    end
    
    %% Save
    save(saveFile,'spikes');
    out = saveFile;
    clear spikes
end
