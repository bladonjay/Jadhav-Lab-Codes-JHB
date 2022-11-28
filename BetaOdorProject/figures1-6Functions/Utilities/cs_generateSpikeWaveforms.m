%cs_getSpikeWaveformsMS
clear
%% Params
% Parameters used in moutainsort
filterRange = [300 6000];
threshold = 3;

% parameters for waveform creation
nPts = 40; % number of points to have in spike waveform
% TODO: Change if mountainsort spike time isn't spike peak
waveWindow = (1:nPts)-14; % puts the actual spike time at pt 14 similar to spike binary where peak is around 14

animals = {'CS39'};
topDir = cs_setPaths;

for a = 1:length(animals)
    animal = animals{a};
    rawDir = 'Z:\Projects\cs_OdorPlace_Project\RAWDataFolders\'; %location of raw data
    animDir = [topDir, animal, 'Expt\',animal,'_direct\']; %location of processed data
    MSDir = 'Z:\Projects\cs_OdorPlace_Project\MountainSort files\'; %location of all MS files
    msDir = [MSDir, animal, '\MountainSort\'];
    
    dayDirs = dir(msDir);
    dayDirs = {dayDirs(3:end).name};
    
    for d = 1:length(dayDirs)
        dayDir = [msDir,dayDirs{d}];
        daynum = str2num(dayDirs{d}(6:7));
        timeFile = dir([dayDir,filesep,'*.timestamps.mda.prv']);
        tmp = jsondecode(fileread([dayDir,filesep,timeFile.name]));
        timeFile = tmp.original_path;
        
        tetDirs = dir(dayDir);
        dr = [tetDirs.isdir];
        tetDirs = {tetDirs(dr).name};
        tetDirs = tetDirs(3:end);
        
        dayStr = getTwoDigitNumber(daynum);
        filePrefix = dayDirs{find(contains(dayDirs,[dayStr,'_']))'};
        filePrefix = strrep(filePrefix,'.mountain','');
        dateStr = filePrefix(end-7:end);
        
        timeDat = [];
        for t = 1:length(tetDirs)
            tetResDir = [dayDir, filesep, tetDirs{t}];
            
            s = [extractBetween(tetResDir, '.nt','.mountain')];
            tet = str2double(s{1});
            %firing and cluster data:
            fireFile = [tetResDir filesep 'firings.curated.mda'];
            if isfile(fireFile)
                
                fireDat = readmda(fireFile); % Rows are # Channels, timestamp (starting @ 0), cluster #
            else
                fireDat = [];
            end
            if isempty(fireDat)
                continue
            end
            
            
            %original path data:
            rawPrv = [tetResDir filesep 'raw.mda.prv'];
            prvDat = jsondecode(fileread(rawPrv));
            rawFile = prvDat.original_path;
            
            %metrics and tags:
            metFile = [tetResDir filesep 'metrics_tagged.json'];
            
            metDat = jsondecode(fileread(metFile));
            metDat = metDat.clusters;
            
            %get clusters
            tags = {metDat.tags};
            rejects = cellfun(@(x) strcmp('rejected',x),tags,'UniformOutput',false);
            acceptcells = ~cellfun(@any, rejects);
            metDat = metDat(acceptcells);
            clusters = [metDat.label];
            
            [clusters,ic] = sort(clusters,'ascend');
            if isempty(clusters)
                continue
            end
            metDat = metDat(ic);
            
            %         pat = '\w+.nt(?<tet>[0-9]+).mda';
            %         parsed = regexp(rawFile,pat,'names');
            %         tet = str2double(parsed.tet);
            disp(['Getting waves from tet ',num2str(tet),'...'])
            
            %         %find original directory path in file name
            %         i = strfind(rawFile,animal);
            %         origDir = rawFile(1:i-1);
            %
            %         %change it to new raw data path
            %         rawFile = strrep(rawFile, origDir ,rawDir);
            %         rawFile = strrep(rawFile,'/','\');
            %
            %         spikeDir = strrep(fileparts(rawFile),'.mda','.spikes');
            
            datePrefix = [dayStr,'_',dateStr];
            dayFolder = [rawDir, animal, '\', datePrefix,'\'];
            %spikeFile = [spikeDir filesep datePrefix '.spikes_nt' num2str(tet) '.dat'];
            spikeFile = [dayFolder,animal,'_',datePrefix,'.spikes\',animal,'_',datePrefix,'.spikes_nt',num2str(tet) '.dat'];
            if ~exist(spikeFile,'file')
                error('Cannot find spikes file at:\n\t%s',spikeFile)
            end
            
            spikeDat = readTrodesExtractedDataFile(spikeFile);
            
            % Load required data
            rawFile = [dayFolder, animal,'_',datePrefix,'.mda\',animal,'_',datePrefix,'.nt',num2str(tet),'.mda'];
            rawDat = readmda(rawFile);
            
            timeFile = strrep(rawFile,['.nt',num2str(tet),'.mda'],'.timestamps.mda');
            if isempty(timeDat)
                timeDat = readmda(timeFile);
            end
            num_channels = size(rawDat,1);
            
            %get time index for all used spikes
            fireIdx = fireDat(2,:);
            fireTimes = timeDat(fireIdx);
            if isrow(fireTimes)
                fireTimes = fireTimes';
            end
            
            tmpFields(1) = struct('name','time','type','uint32','columns',1,'bytesPerItem',4,'data',uint32(fireTimes));
            
            %get wave data for all spikes
            for k=1:num_channels
                
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
            
            %get average waveform across channels
            wavedata = cat(3, tmpFields(2:end).data);
            wavedata = mean(wavedata,3);
            
            %separate waves into clusters
            for c = 1:length(clusters)
                clust = clusters(c);
                ind = find(fireDat(3,:) == clust);
                clustwaves = wavedata(ind,:);
                waves{daynum}{1}{tet}{c} = clustwaves;
            end
            
        end
        
        %save in animal folder
        dayStr = getTwoDigitNumber(daynum);
        save([animDir,animal,'waves',dayStr],'waves');
        clear waves
    end
end

