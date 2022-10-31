%% METADATA

% have Trodes, usrlocal, Src_matlab, matclust, TrodesToMatlab, and Pipeline in path
% To extract matclust files for clustering, run createAllMatclustFiles.m (located in TrodesToMatlab)

% path and animal metadata
topRawDir= 'D:\OdorPlaceAssociation\ThermocoupleTest\CS46Expt\CS46\'; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir= 'D:\OdorPlaceAssociation\ThermocoupleTest\CS46Expt\CS46_direct\';    % data directory
animID= 'CS46';     % your animals name




% tetrode metadata (which tetrode is in which area, which is reference)
TC = 1; 


 % experiment metadata
 cd(topRawDir);
 rawDir=dir();
    rawDir= {rawDir(3:end).name};
 numDays= length(rawDir);
 % The above code assumes raw directory contains one folder for each day of experiment ONLY

 
 
 %% DAY DEPENDENT NON EEG- pos, lfp, spikes, dio, task, linpos
 for sessionNum= 1:numDays
    
%     
    disp(['Creating LFP day ', num2str(sessionNum)]);
    mcz_createNQLFPFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)

    filtPath= ['D:\DataAnalysis\usrlocal\filtering\'];
    addpath('D:\OdorPlaceAssociation\ThermocoupleTest\CS46Expt\CS46_direct\EEG');

    mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);

 end
 createtetinfostruct(dataDir,animID);

sj_addtetrodedescription(dataDir,animID,TC,'TC');

