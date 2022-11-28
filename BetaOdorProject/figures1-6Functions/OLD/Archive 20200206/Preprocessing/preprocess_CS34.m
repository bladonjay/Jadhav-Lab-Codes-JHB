%% METADATA
% have Trodes, usrlocal, Src_matlab, matclust, TrodesToMatlab in path

animID= 'CS34';    % your animal's name

topDir = cs_setPaths(); %set this manually if not on Claire's computer

filtPath= 'D:\DataAnalysis\usrlocal\filtering\'; %set manually - path to usrlocal\filtering folder


topRawDir= [topDir,animID,'Expt\',animID,'\']; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir = [topDir,animID,'Expt\',animID, '_direct\'];

% tetrode metadata (which tetrode is in which area, which is reference)
 CA1   = [1 3 4 17 18 19 20 28 29 30 31 32];
 hpcRef=  2;
 PFC   = [7 8 9 10 13 14 15 16 21 22 23 24 25 26];
 pfcRef=   11;
 OB = 12;
 obRef = 12;
 
 
 % experiment metadata
 cd(topRawDir);
 rawDir=dir();
    rawDir = rawDir([rawDir.isdir]);
    rawDir= {rawDir(3:end).name};
 numDays= length(rawDir);
 % The above code assumes raw directory contains one folder for each day of experiment ONLY
 

 %% DAY DEPENDENT NON EEG- pos, lfp, spikes, dio, task, linpos
 for sessionNum= 1:4
     cs_addPosToSpikesFiles(animID, dataDir, sessionNum)
%      mcz_createNQPosFiles_new([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, [])
%      mcz_createNQLFPFiles([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%      mcz_createNQSpikesFiles([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum) 
%      mcz_createNQDIOFilesFromStateScriptLogs([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%  
 end
   
 % click 4 times on each figure as they pop up: NP --> L well --> choice
 % point --> Rwell
 %              
%        1
%        |
%        |
%        |
%   2----3----4
 
 createtaskstruct(   dataDir,animID, [1 2; 2 2; 2 4; 3 2; 3 4; 4 2], 'getcoord_Tmaze'); %index - [day1 epoch1; day1 epoch2; day2 epoch 1; etc]
 

 
 for sessionNum = 1:8
%      sj_lineardayprocess(dataDir,animID,sessionNum) %get linpos
%      sj_updatetaskstruct(dataDir,animID,sessionNum,[2,4], 'run'); 
%      sj_updatetaskenv(   dataDir,animID,sessionNum,[2,4], 'odorplace');
%      
%       riptet = [20,29,30];
      sj_extractripples(dataDir, animID, sessionNum, riptet, 0.015, 2)
 end
 

 %% DAY INDEPENDENT- cellinfo, tetinfo

 % FF CELL AND TET METADATA STRUCTURES
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 createtetinfostruct(dataDir,animID);
 mcz_createcellinfostruct(dataDir,animID);
 cs_addNumCells(dataDir,animID);
 
%% CELL AND TET Descriptions
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cs_addriptet(dataDir,animID)
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodedescription(dataDir,animID,obRef,'obRef');


 %% EEG PRIMARY PIPELINE

 %array indicating which tetrode should be used as reference. 0 indicates
 %tetrode was missing, or was a reference. 
 refDataTemplate = [2 0 2 2 0 0 11 11 11 11 0 0 11 11 11 11 2 2 2 2 11 11 11 11 11 11 0 2 2 2 2 2];

refData = {};
for sessionNum = 1:numDays
    dayrefData = repmat(refDataTemplate,epochIndex(sessionNum),1);
    refData{1,sessionNum} = dayrefData;
end

% Create all eeg files- raw, and filtered
for sessionNum=1:numDays

 mcz_createRefEEG([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, refData{sessionNum});
%      % EEG wrt ground only
 disp(['doing day',num2str(sessionNum), 'delta']);
mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'deltafilter.mat']);
% 
disp(['doing day',num2str(sessionNum), 'theta']);
 mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
% 
 disp(['doing day ',num2str(sessionNum), ' beta']);
 mcz_betadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'betafilter.mat']);

%      % EEG wrt ground and local reference
disp(['doing day ',num2str(sessionNum), ' gamma']);
mcz_gammadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'gammafilter.mat']);

disp(['doing day ',num2str(sessionNum), ' ripple']);
mcz_rippledayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'ripplefilter.mat']);

end
% 

%% Mean Spectra for each day - 0-40 Hz spectra most useful for beta
%sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB], 1,'fpass',[0 10],'movingwin',[2000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB], 1,'fpass',[0 40],'movingwin',[1000 20]/1000);
%sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB, hpcRef, pfcRef], 0,'fpass',[0 40],'movingwin',[1000 20]/1000);
%sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB, hpcRef, pfcRef], 0,'fpass',[0 100],'movingwin',[400 20]/1000);

