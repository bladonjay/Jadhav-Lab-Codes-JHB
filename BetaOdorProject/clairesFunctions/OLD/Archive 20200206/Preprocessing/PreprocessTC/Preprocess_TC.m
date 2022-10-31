%% METADATA

% have Trodes, usrlocal, Src_matlab, TrodesToMatlab in path
animID= 'CS45';

topDir = ___ %put your top Data directory here (folder containing ThermocoupleTest)
% path and animal metadata
topRawDir= [topDir, 'ThermocoupleTest\', animal, 'Expt\',animal,'\']; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir= [topDir, 'ThermocoupleTest\', animal, 'Expt\',animal,'_direct\'];    % data directory
     % your animals name
 filtPath= 'D:\DataAnalysis\usrlocal\filtering\'; %CHANGE THIS FOR YOUR COMPUTER

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

   
    addpath([dataDir,'EEG\']);

    mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
    mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'deltafilter.mat']);

 end
 
 %% CREATE TETINFO
 createtetinfostruct(dataDir,animID);

sj_addtetrodedescription(dataDir,animID,TC,'TC');

