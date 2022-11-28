%% METADATA

% have Trodes, Src_matlab, TrodesToMatlab in path

% path and animal metadata
topRawDir= ['D:\OdorPlaceAssociation\ThermocoupleTest\CS45Expt\CS45\']; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir= 'D:\OdorPlaceAssociation\ThermocoupleTest\CS45Expt\CS45_direct\';    % data directory
animID= 'CS45';     % your animals name


% tetrode metadata (which tetrode is in which area)
TC = 1; 


 % experiment metadata
 cd(topRawDir);
 rawDir=dir();
    rawDir= {rawDir(3:end).name};
 numDays= length(rawDir);
 % The above code assumes raw directory contains one folder for each day of experiment ONLY

 
 
 %% DAY DEPENDENT - create eeg files, and filtered lfp files
 for sessionNum= 1:numDays
    
%     
    disp(['Creating LFP day ', num2str(sessionNum)]);
    mcz_createNQLFPFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)

    filtPath= ['D:\DataAnalysis\usrlocal\filtering\'];
    addpath('D:\OdorPlaceAssociation\ThermocoupleTest\CS45Expt\CS45_direct\EEG');

    mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
    mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);

 end
 
 %% CREATE TETINFO
 createtetinfostruct(dataDir,animID);

sj_addtetrodedescription(dataDir,animID,TC,'TC');

