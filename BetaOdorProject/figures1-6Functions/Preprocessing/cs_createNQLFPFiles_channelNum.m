function cs_createNQLFPFiles_channelNum(rawDir, dataDir,animID,sessionNum, tet)
%createNQLFPFiles(rawDir, dataDir,animID,sessionNum)
%
%This function extracts LFP information and saves data in the FF format.
%It is assumed that there is
%one (and only one) folder named '*.LFP' in the current working
%directory that contains binary files for each LFP channel (result of extractLFPBinaryFiles.m)  
%
%The function also assumes that there is a '*.time' folder in the current directory
%conaining time information about the session (from
%extractTimeBinaryFile.m).
%
%
%rawDir -- the directory where the raw dat folders are located
%dataDir -- the directory where the processed files should be saved
%animID -- a string identifying the animal's id (appended to the
%beginning of the files).
%sessionNum -- the session number (in chronological order for the animal)
cd(rawDir);

currDir = pwd;
filesInDir = dir;
targetFolder = [];
for i=3:length(filesInDir)
    if filesInDir(i).isdir && ~isempty(strfind(filesInDir(i).name,'.LFP'))
        targetFolder = filesInDir(i).name;
        break;
    end
end

if isempty(targetFolder)
    error('LFP folder not found in this directory.');
end

epochList = getEpochs(1);  %assumes that there is at least a 1-second gap in data between epochs

cd(targetFolder);
datFiles = dir(['*.LFP_nt',num2str(tet),'*.dat']);

if (isempty(datFiles))
    cd(currDir);
    error('No LFP binary files found in LFP folder.');
end

timeDatFiles = dir('*.timestamps.dat');

if (isempty(datFiles))
    cd(currDir);
    error('No timestamps file found in LFP folder.');
end
timeData = readTrodesExtractedDataFile(timeDatFiles(1).name);
timeData = double(timeData.fields(1).data) / timeData.clockrate;

for datFileInd = 1:length(datFiles)
    disp(datFiles(datFileInd).name);
    data = readTrodesExtractedDataFile(datFiles(datFileInd).name);
    
    nTrodeNum = data.ntrode_id;
    channelNum = data.ntrode_channel;
    nTrodeString = getTwoDigitNumber(nTrodeNum);
    channelString = getTwoDigitNumber(channelNum);
        
    for e = 1:size(epochList,1)
        epochString = getTwoDigitNumber(e);
        currentSession = sessionNum;
        sessionString = getTwoDigitNumber(sessionNum);
        currentTimeRange = epochList(e,:);
        
        eeg = []; % mcz from Binary2FF_LFP.m
        epochDataInd = find((timeData >= currentTimeRange(1))&(timeData < currentTimeRange(2)));
        
        eeg{currentSession}{e}{nTrodeNum}.descript = data.description; % mcz from Binary2FF_LFP.m
        eeg{currentSession}{e}{nTrodeNum}.timerange = [timeData(epochDataInd(1)) timeData(epochDataInd(end))]; 
        %redundant notation
        eeg{currentSession}{e}{nTrodeNum}.clockrate = data.clockrate;  % mcz from Binary2FF_LFP.m
        eeg{currentSession}{e}{nTrodeNum}.starttime = timeData(epochDataInd(1));
        eeg{currentSession}{e}{nTrodeNum}.endtime = timeData(epochDataInd(end));       
        eeg{currentSession}{e}{nTrodeNum}.samprate = data.clockrate/data.decimation;
        eeg{currentSession}{e}{nTrodeNum}.nTrode = nTrodeNum;
        eeg{currentSession}{e}{nTrodeNum}.nTrodeChannel = channelNum;
        eeg{currentSession}{e}{nTrodeNum}.data = double(data.fields(1).data(epochDataInd,1)) * data.voltage_scaling *-1;
        eeg{currentSession}{e}{nTrodeNum}.data_voltage_scaled = 1;  % mcz from Binary2FF_LFP.m
        eeg{currentSession}{e}{nTrodeNum}.data_voltage_inverted = 1;  
        
        if strcmp(data.reference, 'on')                      % mcz from Binary2FF_LFP.m
            eeg{currentSession}{e}{nTrodeNum}.referenced = 1; % mcz from Binary2FF_LFP.m
        else                                                  % mcz from Binary2FF_LFP.m
            eeg{currentSession}{e}{nTrodeNum}.referenced = 0; % mcz from Binary2FF_LFP.m
        end                                                   % mcz from Binary2FF_LFP.m
        
        eeg{currentSession}{e}{nTrodeNum}.low_pass_filter = data.low_pass_filter; % mcz from Binary2FF_LFP.m
        eeg{currentSession}{e}{nTrodeNum}.voltage_scaling = data.voltage_scaling; % mcz from Binary2FF_LFP.m


        
        lfpdir = pwd;
        cd(dataDir);
        %mkdir EEG
        save([dataDir,'EEG/',animID,'eeg',sessionString,'-',epochString,'-',nTrodeString,'ch',channelString,'.mat'],'-v6','eeg');
        cd(lfpdir);
    end
end
cd(currDir);    

function numString = getTwoDigitNumber(input)
    
if (input < 10)
    numString = ['0',num2str(input)];
else
    numString = num2str(input);
end
    
