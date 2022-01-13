function [units,out] = convert_ml_to_struct(dayDir,varargin)
%animDir is animal_direct folder with processed data.

%sessionNum is day number

samplerate = 30000;
min_epoch_break=4; % this should be flexible, is currently unused
overwrite = 1;

assignVars(varargin)

%% Find correct day folder inside MountainSort results
% need to find where the 'MountainSort directory is:
if contains(dayDir,'.mountain')
    [~,sessionstr]=fileparts(dayDir);
    pat='(?(?<anim>[A-Z]+[0-9]+)_(?<sesn>[^_]+)_(?<date>[^_]+)';
    sessInfo=regexp(sessionstr,pat,'names');
    daystr=[sessInfo.sesn '_' sessInfo.date];
else
    fprintf('Need to find a folder with the ending .mountain that contains tetrode folders!')
    return
end




%% Collect Data and create save file
tetDirs = dir([dayDir filesep '*nt*']);
timeFile = dir([dayDir filesep '*timestamps*']);
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
        fprintf('I couldnt find time file, can you?\n ');
        timeDir=uigetdir(dayDir,'Find the .timestamps file');
        timeFile = dir([timeDir filesep '*timestamps*']);
        if ~isempty(timeFile)
            timeDat = readmda(fullfile(timeFile.folder,timeFile.name));
        else
            fprintf('You couldnt find the file either, quitting!\n');
            return
        end
    end
end
% need to update file naming convention so that its not just sessnum
% but also date!!!

saveFile = sprintf('%s_spikes.mat',sessionstr);
saveDir=fileparts(dayDir);
fprintf('saving file %s \n    in dir %s \n',saveFile, saveDir);
if exist(fullfile(saveDir,saveFile),'file') && ~overwrite
    fprintf('Oops!!! \n     Spikes file @ %s already exists. Skipping...\n',fullfile(saveDir,saveFile));
    return;
end

%% get tet numbers
pat = '\w*.nt(?<tet>[0-9]+).\w*';
parsed = cellfun(@(x) regexp(x,pat,'names'),{tetDirs.name});
tet_nums = str2double({parsed.tet});
Ntets = max(tet_nums);
fprintf('Detected %i tetrodes. Tetrodes:\n',Ntets);
disp(tet_nums);
drawnow;

% this is in essence filter framework
%spikes = cell(1,sessionNum);
%spikes{sessionNum} = cell(1,Nepochs);
%spikes{sessionNum}(:) = {cell(1,Ntets)};

% this is the struct composition:
% units is struct with fields:
% ts
% sessionNum
% tetrode
% unitname (a-z on that tet)
% etc
units=struct();

it=1; % we dont know how many cells there are yet, so we'll iterate!
%% Convert from firings data
for k=1:numel(tetDirs)
    tD = [tetDirs(k).folder filesep tetDirs(k).name];
    tetNum = tet_nums(k);
    metFile = [tD filesep 'metrics_tagged.json'];
    
    %tetrodes with no cells will not have firings.curated.mda files
    if isfile([tD filesep 'firings.curated.mda'])
        fireFile = [tD filesep 'firings.curated.mda'];
    else
        fireFile = [tD filesep 'firings_raw.mda'];
    end
    metDat = jsondecode(fileread(metFile));
    metDat = metDat.clusters;
    
    %CS- remove cells that are tagged as "rejected." Keep all other tags
    %and save in spikes file.
    tags = {metDat.tags};
    rejects = cellfun(@(x) strcmp('rejected',x),tags,'UniformOutput',false);
    acceptcells = ~cellfun(@any, rejects);
    metDat = metDat(acceptcells);
    
    fireDat = readmda(fireFile); % Rows are # Channels, timestamp (starting @ 0), cluster #
    
    clusters = [metDat.label];
    %Nclust = max(clusters)+any(clusters==0);
    [clusters,ic] = sort(clusters,'ascend');
    metDat = metDat(ic);
    
    % TODO: check params.json for smaplerate and error if not matching
    
    fireTimes = timeDat(fireDat(2,:));
    unitID='A';
    for m=1:length(clusters)
        mySpikes=fireTimes(fireDat(3,:)==clusters(m))';
        mD = metDat(m).metrics;
        units(it).fields = 'spikeInd spikeTime MSepochDetected';
        units(it).info = 'spike data from MountainSort 4 (MountainLab-JS)';
        units(it).fs =samplerate;
        units(it).meanrate = mD.firing_rate;
        units(it).isolation = mD.isolation;
        units(it).num_events = mD.num_events;
        units(it).noise_overlap = mD.noise_overlap;
        units(it).peak_noise = mD.peak_noise;
        units(it).peak_amplitude = mD.peak_amp;
        units(it).peak_snr = mD.peak_snr;
        units(it).spikes = [mySpikes mySpikes/samplerate fireDat(1,fireDat(3,:)==clusters(m))' ]; % 1
        units(it).tag=metDat(m).tags;
        units(it).tetNum=tetNum;
        units(it).unitID=unitID;
        %% add tag
        
        tag = metDat(m).tags;
        if length(tag) > 1
            tag = tag{(~(cellfun('isempty',strfind(tag,'accepted'))))};
        end
        
        % now iterate up
        it=it+1;
        if double(unitID)==90
            unitID='a'; % skip to lower case if we get to Z
        else
            unitID=char(unitID+1);
        end
        
    end % of units
end % of tetrodes


%% Save
save(fullfile(saveDir,saveFile),'units');
out = fullfile(saveDir,saveFile);
end
