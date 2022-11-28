function out = convert_ml_to_FF(animID,animDir,sessionNum)
%animDir is animal_direct folder with processed data.
%sessionNum is day number

samplerate = 30000;
min_epoch_break=1;
overwrite = 1;

%assignVars(varargin)

%% Find correct day folder inside MountainSort results
dayfolders = dir([animDir filesep 'MountainSort']);
%dayfolders = dir('G:\Data\OdorPlaceAssociation\CS44Expt\CS44_direct\MountainSort');
dayfolders = {dayfolders(3:end).name};

sessionstr = sprintf('%s_%02i',animID,sessionNum);
%resDir = ['G:\Data\OdorPlaceAssociation\CS44Expt\CS44_direct\MountainSort', filesep, dayfolders{contains(dayfolders,sessionstr)}];
resDir = [animDir, filesep, 'MountainSort', filesep, dayfolders{contains(dayfolders,sessionstr)}];


%% Collect Data and create save file
out = '';
tetDirs = dir([resDir filesep '*.mountain']);
timeFile = dir([resDir filesep '*timestamps.mda']);
if isempty(timeFile)
    error('No timestamps file found in %s',resDir);
end
timeFile = [timeFile.folder filesep timeFile.name];
if strcmp(timeFile(end-2:end),'prv')
    tmp = jsondecode(fileread(timeFile));
    timeFile = tmp.original_path;
end
timeDat = readmda(timeFile);
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
    if isfile([tD filesep 'firings.curated.mda'])
        fireFile = [tD filesep 'firings.curated.mda'];
        %fireFile2 = [tD filesep 'firings_raw.mda'];
        metDat = jsondecode(fileread(metFile));
        metDat = metDat.clusters;
        
        %CS- remove cells that are tagged as "rejected." Keep all other tags
        %and save in spikes file.
        tags = {metDat.tags};
        rejects = cellfun(@(x) strcmp('rejected',x),tags,'UniformOutput',false);
        acceptcells = ~cellfun(@any, rejects);
        metDat = metDat(acceptcells);
        
        fireDat = readmda(fireFile); % Rows are # Channels, timestamp (starting @ 0), cluster #
        %fireDat2 = readmda(fireFile2);
        
        clusters = [metDat.label];
        %Nclust = max(clusters)+any(clusters==0);
        [clusters,ic] = sort(clusters,'ascend');
        if isempty(clusters)
            continue
        else
        metDat = metDat(ic);
        
        % TODO: check params.json for smaplerate and error if not matching
        
        fireTimes = timeDat(fireDat(2,:));
        for ep=1:Nepochs
            spikes{sessionNum}{ep}{tetNum} = cell(1,length(clusters));
            for c=1:length(clusters)
                t1 = epoch_starts(ep);
                if ep<Nepochs
                    t2 = epoch_starts(ep+1)-1;
                else
                    t2=timeDat(end);
                end
                
                
                idx = fireDat(3,:)==clusters(c) & fireTimes>=t1 & fireTimes<=t2;
                mD = metDat(c).metrics;
                dat = zeros(sum(idx),7);
                fields = 'time x y dir not_used amplitude(highest variance channel) posindex n_detection_channels';
                descript = 'spike data from MountainSort 4 (MountainLab-JS)';
                meanrate = mD.firing_rate;
                peak_amplitude = mD.peak_amp;
                %refractory_violation_1msec = mD.refractory_violation_1msec;
                dat(:,1) = fireTimes(idx)'/samplerate;
                dat(:,7) = fireDat(1,idx)';
                %% add tag
                
                if isempty(dat)
                    %spikes{sessionNum}{ep}{tetNum}{c} = [];
                else
                    
                    tag = metDat(c).tags;
                    if length(tag) > 1
                        tag = tag{(~(cellfun('isempty',strfind(tag,'accepted'))))};
                    end
                    
                    spikes{sessionNum}{ep}{tetNum}{c} = struct('data',dat,'meanrate',meanrate,...
                        'descript',descript,'fields',fields,...
                        'timerange',[t1 t2]/samplerate,...
                        'tag',tag,...
                        'peak_amplitude',peak_amplitude);
                end
                
            end
        end
        end
    end
end

%% Save
save(saveFile,'spikes');
out = saveFile;
