function out = md_convert_ml_to_FF(animID,resDir, saveDir, sessionNum,varargin)
% Convert mountainsort curated spikes to filter framework spikes files
% animID : animal ID
% resDir : mountainlab results directory for a single day
% sessionNum: number of day
% NAME-VALUE Pairs:
% samplerate : sampling rate in Hz, default 30000
% min_epoch_break : minumum recording time break to classify an epoch, 1 sec default
% overwrite : 1 or 0 (default), whether to overwrite an existing spikes file

% saves spikes file. Only times (column 1) and # channels (col 7)

samplerate = 30000;
min_epoch_break=1;
overwrite = 0;

assignVars(varargin)

out = '';
tetDirs = dir([resDir filesep '*.mountain']);
timeFile = dir([resDir filesep '*timestamps*']);
posFile = dir([saveDir filesep animID 'pos' sprintf('%02i', sessionNum) '.mat']);

if isempty(timeFile)
    error('No timestamps file found in %s',resDir);
end
timeFile = [timeFile.folder filesep timeFile.name];
if strcmp(timeFile(end-2:end),'prv')
    tmp = jsondecode(fileread(timeFile));
    timeFile = tmp.original_path;
end
timeDat = readmda(timeFile);
saveFile = sprintf('%s%s%sspikes%02i.mat',saveDir,filesep,animID,sessionNum);
if exist(saveFile,'file') && ~overwrite
    fprintf('Spikes file @ %s already exists. Skipping...\n',saveFile)
    return;
end

% Determine epoch start times by looking for gaps longer than 1 sec
gaps = diff(timeDat);
epoch_gaps = find(gaps>=min_epoch_break*samplerate);
epoch_starts = timeDat([1 epoch_gaps])./samplerate; % in sec
Nepochs = numel(epoch_starts);
fprintf('Detected %i epochs. Start times:\n',Nepochs);
disp(epoch_starts');

track_start = zeros(Nepochs, 1); % time of tracking started in each epoch
track_end = Inf(Nepochs, 1); % time of tracking ended
if ~isempty(posFile)
    load([posFile.folder filesep posFile.name], 'pos');
    for epoch = 1:length(pos{sessionNum})
        track_start(epoch) = pos{sessionNum}{epoch}.data(1,1);
        track_end(epoch) = pos{sessionNum}{epoch}.data(end,1);
    end
else
    warning('No position file found in %s', saveDir);
end

% get tet numbers
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

for k=1:numel(tetDirs)
    tD = [tetDirs(k).folder filesep tetDirs(k).name];
    tetNum = tet_nums(k);
    metFile = [tD filesep 'metrics_tagged.json'];
    fireFile = [tD filesep 'firings.curated.mda'];
    if ~isfile(fireFile)
        fprintf('No curated firings file found for tetrode %i. Skipping...\n',tetNum);
        continue;
    end
    
    metDat = jsondecode(fileread(metFile));
    metDat = metDat.clusters;
    fireDat = readmda(fireFile); % Rows are # Channels, timestamp (starting @ 0), cluster #
    
    clusters = [metDat.label];
    %         Nclust = max(clusters)+any(clusters==0);
    Nclust = length(clusters);
    
    [clusters,ic] = sort(clusters,'ascend');
    metDat = metDat(ic);
    
    % TODO: check params.json for samplerate and error if not matching
    
    fireTimes = timeDat(fireDat(2,:)+1)./samplerate; % convert to time in sec % +1 given timestamp starts at 0
    for l=1:Nepochs
        spikes{sessionNum}{l}{tetNum} = cell(0);
        t1 = max([epoch_starts(l) track_start(l)]);
        if l<Nepochs
            t2 = min([epoch_starts(l+1)-1 track_end(l)]);
        else
            t2= min([timeDat(end)/samplerate track_end(l)]);
        end
        for m=1:Nclust
            tag = metDat(m).tags;
            if any(strcmp(tag,'accepted'))
                tag = tag{strcmp(tag,'accepted')};
                %                 if iscell(tag)
                %                     tag = cell2mat(tag);
                idx = fireDat(3,:)==clusters(m) & fireTimes>=t1 & fireTimes<=t2;
                mD = metDat(m).metrics;
                dat = zeros(sum(idx),8);
                fields = 'time x y dir not_used amplitude(highest variance channel) posindex n_detection_channels';
                descript = 'spike data from MountainSort 4 (MountainLab-JS)';
                meanrate = mD.firing_rate;
                peak_amplitude = mD.peak_amp;
                % refractory_violation_1msec = mD.refractory_violation_1msec; % No longer found in curated mda, possibly due to version update
                if sum(idx)
                    dat(:,1) = fireTimes(idx)';
                    posindex = lookup(dat(:,1), pos{sessionNum}{l}.data(:,1));
                    dat(:,2:4) = pos{sessionNum}{l}.data(posindex, 2:4);
                    dat(:,7) = posindex;
                    %                 if
                    %                 data(:,2:4) =
                    
                    dat(:,8) = fireDat(1,idx)';
                    %             tag = metDat(m).tags;
                    %             if any(strcmp(tag,'accepted'))
                    %                 tag = tag{strcmp(tag,'accepted')};
                    %             end
                else
                    dat = [];
                end
                spikes{sessionNum}{l}{tetNum}{end+1} = struct('data',dat,'meanrate',meanrate,...
                    'descript',descript,'fields',fields,...
                    'timerange',[t1 t2],...
                    'tag',tag,...
                    'peak_amplitude',peak_amplitude,...
                    'clusterID', clusters(m));
            end
        end
    end
end

save(saveFile,'spikes');
out = saveFile;
