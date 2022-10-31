function out = cs_generateSpikeWaveforms(topDir, animal, day, tet)
%CS- modified from 'generateMatclustFromMountainSort,' this function will
%add spike waveforms to spikes files after sorting and curation in
%MountainSort. Requires extracted spikes.dat files as well as MountainSort files.

% Parameters used in moutainsort
filterRange = [300 6000];
threshold = 3;

% parameters for waveform creation
nPts = 40; % number of points to have in spike waveform
% TODO: Change if mountainsort spike time isn't spike peak
waveWindow = (1:nPts)-14; % puts the actual spike time at pt 14 similar to spike binary where peak is around 14

rawDir = [topDir, animal,'Expt\',animal,'\'];
%dataDir = [topDir, animal,'Expt\',animal,'_direct\'];
dataDir = 'G:\Data\OdorPlaceAssociation\CS44Expt\CS44_direct\';

daystr = getTwoDigitNumber(day);
dayname = dir([rawDir,daystr,'*']);
dayname = dayname(1).name;

tetResDir = [dataDir,'MountainSort\',animal,'_',dayname,'.mountain\',animal,'_',dayname,'.nt',num2str(tet),'.mountain'];

% Get all relevant file names
fireFile = [tetResDir filesep 'firings.curated.mda'];
fireFile = [tetResDir filesep 'firings_raw.mda'];

fprintf('Loading mountainsort data...\n')
% Get and load trodes extracted spikes binary


spikeDir = dir([rawDir,dayname,'\',animal,'*.spikes']);
spikeDir = [spikeDir(1).folder, filesep, spikeDir(1).name];
[~,datPrefix] = fileparts(spikeDir);

spikeFile = [spikeDir filesep datPrefix '.spikes_nt' num2str(tet) '.dat'];
if ~exist(spikeFile,'file')
    error('Cannot find spikes file at:\n\t%s',spikeFile)
end
spikeDat = readTrodesExtractedDataFile(spikeFile);

mdaDir = dir([rawDir,dayname,'\',animal,'*.mda']);
mdaDir = [mdaDir(1).folder, filesep, mdaDir(1).name];
%rawFile = strrep(rawFile,'/media/roshan/ExtraDrive1','/data');
rawFile = [mdaDir,filesep,datPrefix,'.nt',num2str(tet),'.mda'];

timeFile = [mdaDir, filesep,datPrefix,'.timestamps.mda'];
metricsFile = [tetResDir filesep 'metrics_tagged.json'];

fireFile = [tetResDir filesep 'firings.curated.mda'];

% Load required data
rawDat = readmda(rawFile);
timeDat = readmda(timeFile);
fireDat = readmda(fireFile);
metricsDat = jsondecode(fileread(metricsFile));

clusttags = {metricsDat.clusters.tags};
acceptcells = ~cellfun(@any, cellfun(@(x) strcmp('rejected',x),clusttags,'UniformOutput',false));
clusters = [metricsDat.clusters(acceptcells).label];

test = fireDat(:,ismember(fireDat(3,:),clusters));

spikeStruct.num_channels = size(rawDat,1);
%Add code to determine if there are dead channels- exclude these



% create fields
fprintf('Gathering spike waveforms...\n')
fireIdx = fireDat(2,:);
fireTimes = timeDat(fireIdx);
if isrow(fireTimes)
    fireTimes = fireTimes';
end
tmpFields(1) = struct('name','time','type','uint32','columns',1,'bytesPerItem',4,'data',uint32(fireTimes));
for k=1:spikeStruct.num_channels
    tmpFields(k+1).name = ['waveformCh' num2str(k)];
    tmpFields(k+1).type = 'int16';
    tmpFields(k+1).columns = nPts;
    waveDat = int16(zeros(numel(fireIdx),nPts));
    for l=1:numel(fireIdx)
        tmpWin = fireIdx(l)+waveWindow;
        
        % Fix windows at start of data
        if any(tmpWin<=0)
            tmpWin = 1:nPts;
        end
        
        % Fix windows at end of data
        if any(tmpWin>size(rawDat,2))
            tmpWin = size(rawDat,2)-nPts+1:size(rawDat,2);
        end
        
        waveDat(l,:) = int16(rawDat(k,tmpWin));
    end
    tmpFields(k+1).data = waveDat;
end
spikeStruct.fields = tmpFields;
data = spikeStruct;
fprintf('Waveforms gathered.\n')
clearvars -except data dataDir tetResDir fireDat timeDat metricsDat



%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
%% data structure creation complete. Now can use code from
%  createMatclustFile to generate matclust param and waves files
%--------------------------------------------------------------------------------
%--------------------------------------------------------------------------------
fprintf('Generating matclust param and waves files...\n')
maxSamplesToPeak = 17;
minSamplesToPeak = 10;

spikeTimes = double(data.fields(1).data)/data.clockrate;
spikeClusters = fireDat(3,:);

%put all spikes into one 3d matrix ([spikeind spikelength channel])
spikeWaveforms = zeros(length(spikeTimes),data.fields(2).columns,length(data.fields)-1);
peakVals = [];
peakInds = [];
for fieldInd = 2:length(data.fields)
    spikeWaveforms(:,:,fieldInd-1) = double(data.fields(fieldInd).data)*data.voltage_scaling;
    [tmpPeakVals, tmpPeakInds] = max(spikeWaveforms(:,:,fieldInd-1),[],2);
    peakVals = [peakVals tmpPeakVals];
    peakInds = [peakInds tmpPeakInds];
end

[junk, maxChannels] = max(peakVals,[],2);

keeperPeakInds = zeros(length(maxChannels),1);
for i=1:length(maxChannels)
    keeperPeakInds(i) = peakInds(i,maxChannels(i));
end
keeperPeakInds = find((keeperPeakInds <= maxSamplesToPeak)&(keeperPeakInds >= minSamplesToPeak));
spikeWaveforms = spikeWaveforms(keeperPeakInds,:,:);
spikeClusters = spikeClusters(keeperPeakInds);
spikeTimes = spikeTimes(keeperPeakInds);

windowSize = size(spikeWaveforms,2);
midpoint = round(windowSize/2);

if (size(spikeWaveforms,1) > 100)
    
    scores = [];
    haveStatsToolbox = which('princomp');
    
    %if the user has access to the 'pca' function, calculate the 1st and
    %2nd priciple components
    %         if ~isempty(haveStatsToolbox)
    %             for ch = 1:size(spikeWaveforms,3)
    %                 pcaWaves = spikeWaveforms(:,:,ch);
    %                 for i = 1:size(pcaWaves,2)
    %                     pcaWaves(:,i) = pcaWaves(:,i)-mean(pcaWaves(:,i));
    %                 end
    %
    %                 %coef = pca(pcaWaves);
    %                 coef = princomp(pcaWaves);
    %                 coef = coef(:,1:2); %first two components
    %                 %coef = coef(:,1); %just keep the 1st component
    %                 scores = [scores (pcaWaves*coef)];
    %
    %             end
    %         end
    filedata.paramnames = {'Time'};
    filedata.params = [spikeTimes];
    for i = 1:size(spikeWaveforms,3)
        filedata.paramnames = [filedata.paramnames,['Peak ',num2str(i),' (uV)']];
        filedata.params = [filedata.params max(spikeWaveforms(:,:,i),[],2)];
    end
    for i = 1:size(spikeWaveforms,3)
        filedata.paramnames = [filedata.paramnames,['Trough ',num2str(i),' (uV)']];
        filedata.params = [filedata.params min(spikeWaveforms(:,:,i),[],2)];
    end
    for i = 1:size(spikeWaveforms,3)
        filedata.paramnames = [filedata.paramnames,['Peak to trough ',num2str(i),' (uV)']];
        filedata.params = [filedata.params max(spikeWaveforms(:,:,i),[],2)-min(spikeWaveforms(:,:,i),[],2)];
    end
    
    if ~isempty(haveStatsToolbox)
        for i = 1:size(spikeWaveforms,3)
            filedata.paramnames = [filedata.paramnames,['1stPCA ',num2str(i)]];
            filedata.paramnames = [filedata.paramnames,['2ndPCA ',num2str(i)]];
        end
    end
    filedata.params = [filedata.params scores];
    
    if isrow(spikeClusters)
        spikeClusters = spikeClusters';
    end
    filedata.params = [filedata.params spikeClusters];
    filedata.paramnames = [filedata.paramnames,'Cluster'];
    
    
    waves = permute(spikeWaveforms,[2 3 1]);
else
    disp('Not enough spike events in file');
    filedata = [];
    waves = [];
    
end

out = [];
%-----------------------------------
if (~isempty(filedata))
    saveName = fullfile(destFolder,['param_nt',num2str(data.ntrode_id)]);
    waveName = fullfile(destFolder,['waves_nt',num2str(data.ntrode_id)]);
    filedata.filename = ['waves_nt' num2str(data.ntrode_id)];
    save(saveName, '-v7.3', 'filedata');
    save(waveName, '-v7.3', 'waves');
    out = {saveName;waveName};
else
    error('No param data')
end

%    %--------------------------------------------------------------------------------
%    %--------------------------------------------------------------------------------
%
%    %% Now hack-and-slash a matclust save file in order to import the moutainsort clustering to matclust
%
%    %--------------------------------------------------------------------------------
%    %--------------------------------------------------------------------------------
%
%    fprintf('Now creating a matclust save file containing the mountainsort cluster...\n')
%    min_epoch_gap = 1;
%    samplerate = data.clockrate;
%
%    % Load matclust save file template
%    mcTemplate = load('matclust_file_template.mat');
%    clustdata = mcTemplate.clustdata;
%    clustattrib = mcTemplate.clustattrib;
%    graphattrib = mcTemplate.graphattrib;
%
%    %% Customize clustdata
%    nParams = size(filedata.params,2);
%    clustdata.filledparam = ones(1,nParams);
%    clustdata.params = filedata.params;
%    clustdata.origparams = nParams;
%    clustdata.names = filedata.paramnames;
%    nSpikes = size(filedata.params,1);
%    spikeTimes = filedata.params(:,1);
%
%    % create time filters
%    if isrow(timeDat)
%        timeDat = timeDat';
%    end
%    epoch_gaps = find(diff(timeDat)>=min_epoch_gap*samplerate);
%    epoch_start = double(timeDat([1 epoch_gaps+1]))./samplerate;
%    epoch_end = double(timeDat([epoch_gaps numel(timeDat)]))./samplerate;
%    timefilterranges = [epoch_start epoch_end];
%    timefilterranges = [double(timeDat([1 end]))'./samplerate+[-550 550];timefilterranges];
%    clustdata.timefilterranges = timefilterranges;
%
%    timefilters = int32(ones(nSpikes,1));
%    clustdata.otherfilters = timefilters;
%    timefilternames = clustdata.timefilternames;
%    for k=2:size(timefilterranges,1)
%        tfr = timefilterranges(k,:);
%        spikes_in_range = spikeTimes>=tfr(1) & spikeTimes<=tfr(2);
%        timefilters(spikes_in_range) = k*2-1;
%        tSec = fix(rem(tfr,60));
%        tMin = fix(rem(tfr/60,60));
%        tHour = fix(tfr/60/60);
%        timefilternames{k} = sprintf('%i e%i %i:%02i:%02i-%i:%02i:%02i',k,k-1,tHour(1),tMin(1),tSec(1),tHour(2),tMin(2),tSec(2));
%    end
%    clustdata.timefilters = timefilters;
%
%    clustdata.timefiltermemmap = zeros(32,1);
%    clustdata.timefiltermemmap(1:size(timefilterranges,1)) = 1:size(timefilterranges,1);
%    clustdata.filtermemmap = [clustdata.timefiltermemmap;clustdata.otherfiltermemmap];
%    clustdata.filteredpoints = true(nSpikes,1);
%    clustdata.dataranges = [min(filedata.params);max(filedata.params)];
%
%    %% Customize graphattrib
%    graphattrib.nonhiddenpoints = true(nSpikes,1);
%    graphattrib.currentdatarange = clustdata.dataranges(:,1:2);
%
%    %% Customize clustattrib
%    nClust = numel(unique(spikeClusters));
%    clustNums = sort(unique(spikeClusters),'ascend');
%    clusters = cell(1,nClust);
%    for k=1:nClust
%        clusters{k}.defineaxes=[];
%        clusters{k}.box={};
%        clusters{k}.polyindex=[];
%        clusters{k}.index=find(spikeClusters==clustNums(k));
%    end
%    clustattrib.clusters = clusters;
%    clustattrib.filterindex = {[1 1]};
%    clustattrib.clustersOn = ones(nClust,1);
%    clustattrib.currentfilepath = destFolder;
%    clustattrib.currentfilename = ['mountain_param_nt' num2str(data.ntrode_id) '_template.mat'];
%    clustattrib.currentparamfilename = ['param_nt' num2str(data.ntrode_id) '.mat'];
%    clustattrib.dataFile = ['waves_nt' num2str(data.ntrode_id) '.mat'];
%    clustattrib.pointexclude = int32(zeros(nSpikes,1));
%    clustattrib.pointinclude = int32(zeros(nSpikes,1));
%    clustattrib.currclust = 1;
%
%    save([destFolder filesep clustattrib.currentfilename],'clustdata','graphattrib','clustattrib')
fprintf('Done!\n')







