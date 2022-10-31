%% METADATA

% have Trodes, usrlocal, Src_matlab, matclust, TrodesToMatlab, and Pipeline in path
% To extract matclust files for clustering, run createAllMatclustFiles.m (located in TrodesToMatlab)

% path and animal metadata
animID= 'CS31'; 
topDir = cs_setPaths();
topRawDir= [topDir,animID,'Expt\',animID,'\']; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory
    % your animals name

% tetrode metadata (which tetrode is in which area, which is reference)
 CA1   = [3 4 17 18 19 20 27 28 29 30 31 32];
 hpcRef=  1;
 PFC   = [7 8 9 12 13 14 15 16 21 22 23 25 26];
 pfcRef=   10;
 OB = 24;
 obRef = 24; 
 
refDataTemplate = [0 1 1 1 0 0 10 10 10 0 10 10 10 10 10 10 1 1 1 1 10 10 10 0 10 10 1 1 1 1 1 1];

 
 % experiment metadata
 cd(topRawDir);
 rawDir=dir();
    rawDir= {rawDir(3:end).name};
 numDays= length(rawDir);
 % The above code assumes raw directory contains one folder for each day of experiment ONLY

 %% DAY DEPENDENT NON EEG- pos, lfp, spikes, dio, task, linpos
 for sessionNum= 1:3 %1:numDays;
     
 % CONVERT RAW DATA TO BASIC FF DATA
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
     %mcz_createNQPosFiles_new([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, [])
%   
%      mcz_createNQLFPFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%      mcz_createNQSpikesFiles([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%      mcz_createNQDIOFilesFromStateScriptLogs([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
      cs_addPosToSpikesFiles(animID, dataDir, sessionNum)
 end
 
     % FF BEHAVIOR METADATA STRUCTURES
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      createtaskstruct(   dataDir,animID, [1 2; 1 4; 2 2; 2 4; 2 6; 3 2; 3 4; 3 6], 'getcoord_Tmaze'); %index - [day1 epoch1; day1 epoch2; etc]
     
     % FF SECONDARY DATA STRUCTURES
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for sessionNum= 1:3    
%      sj_lineardayprocess(dataDir,animID,sessionNum) %get linpos
     %sj_updatetaskstruct(dataDir,animID,sessionNum,[#], 'sleep'); % one line for each day
     %sj_updatetaskenv(   dataDir,animID,sessionNum,[2 4 ], '');
     
end

 %% DAY INDEPENDENT- cellinfo, tetinfo

 % FF CELL AND TET METADATA STRUCTURES
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createtetinfostruct(dataDir,animID);
 mcz_createcellinfostruct(dataDir,animID);
 cs_addNumCells(dataDir,animID);

 
 % CELL AND TET DESCRIPTIONS
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %riptet = [17,18,30,32];
%sj_addtetrodedescription(dataDir,animID,riptet,'riptet');

cs_riptettag({animID});
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');
sj_addtetrodedescription(dataDir,animID,obRef,'obRef');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodelocation(dataDir,animID,OB,'OB');
%   

for sessionNum = 1:3
    riptet = cs_getriptet_mostcells(animID,day);
    sj_extractripples(dataDir, animID, sessionNum, CA1, 0.015, 2)
end
 
%  betatets = [7 17 18 21 24 30 32];
%  cssj_addtetrodedescription2(dataDir,'CS31',betatets,'betatet');

 %% EEG PRIMARY PIPELINE

%Creates refData array with different numbers of copies of refDataTemplate
%on different days according to number of epochs- siplifies so script
%doesn't have to be run multiple times for different days. 

%Create epoch index (number of epochs per day)
 epochIndex = [];
 for sessionNum= numDays:-1:1
     cd(topRawDir);
     dayDir = rawDir{sessionNum};
     cd(dayDir)
     recfiles = dir('*.rec');
     numEpochs = length(recfiles);
     epochIndex(sessionNum) = numEpochs;
 end 
 
 
refData = {};
for sessionNum = numDays:-1:1
    dayrefData = repmat(refDataTemplate,epochIndex(sessionNum),1);
    refData{1,sessionNum} = dayrefData;
end
 
 
filtPath= 'D:\DataAnalysis\usrlocal\filtering\';
addpath('D:\OdorPlaceAssociation\CS31Expt\CS31_direct\EEG');
for sessionNum=1:3
    
%mcz_createRefEEG([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, refData);
     
% EEG wrt ground only
%mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'deltafilter.mat']);
mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
disp(['doing day ',num2str(sessionNum), ' beta']);
%mcz_betadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'betafilter.mat']);
    

% EEG wrt ground and local reference
%mcz_gammadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'gammafilter.mat']);
%mcz_rippledayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'ripplefilter.mat']);
end
% 



%sjcs_baselinespecgram(prefix, days, epochs, tetrodelist, do_wrtgnd,'fpass',fpass,'movingwin',movingwin);

% sjcs_baselinespecgram(animID, 1:numDays, [2 4 6], [CA1, PFC, OB], 1,'fpass',[0 10],'movingwin',[2000 20]/1000);
% sjcs_baselinespecgram(animID, 1:numDays, [2 4 6], [CA1, PFC, OB], 1,'fpass',[0 40],'movingwin',[1000 20]/1000);
% sjcs_baselinespecgram(animID, 1:numDays, [2 4 6], [CA1, PFC, OB, hpcRef, pfcRef], 0,'fpass',[0 40],'movingwin',[1000 20]/1000);
% sjcs_baselinespecgram(animID, 1:numDays, [2 4 6], [CA1, PFC, OB, hpcRef, pfcRef], 0,'fpass',[0 100],'movingwin',[400 20]/1000);
