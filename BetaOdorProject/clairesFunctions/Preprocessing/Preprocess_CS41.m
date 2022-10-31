%% METADATA

% have Trodes, usrlocal, Src_matlab, TrodesToMatlab in path
animID= 'CS41';    % your animal's name

topDir = cs_setPaths(); %set this manually if not on Claire's computer

filtPath= 'D:\DataAnalysis\usrlocal\filtering\'; %set manually - path to usrlocal\filtering folder

topRawDir= [topDir,animID,'Expt\',animID,'\']; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir = [topDir,animID,'Expt\',animID, '_direct\'];

% tetrode metadata (which tetrode is in which area)
TC = 24; 
OB = 12;
CA1 = [17 18 19 20 28 29 30 32 1 2 3 4];
hpcRef = 27;
PFC = [5 7 8 10 21 22 23 25 26 11 13 14 15 16];
pfcRef = 9;


 % experiment metadata
%  cd(topRawDir);
%  rawDir=dir();
%     rawDir= {rawDir(3:end).name};
%  numDays= length(rawDir);
numDays = 10;
 % The above code assumes raw directory contains one folder for each day of experiment ONLY

%Create refDataTemplate: array indicating which tetrode should be used as
%reference. 0 indicates tetrode was missing, or was a reference. 
allTets = [TC,OB,CA1,hpcRef,PFC,pfcRef];
refDataTemplate = 1:32;
unused = find(~ismember(refDataTemplate,allTets));
refDataTemplate([TC, OB, hpcRef, pfcRef,unused]) = 0;
refDataTemplate(CA1) = hpcRef;
refDataTemplate(PFC) = pfcRef;
 
 
 %% DAY DEPENDENT - create eeg files, and filtered lfp files
 for sessionNum= 2:numDays
    
%      disp(['Creating pos day ', num2str(sessionNum)]); 
%       mcz_createNQPosFiles_new([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, [])
% %     
%     disp(['Creating LFP day ', num2str(sessionNum)]);
%     mcz_createNQLFPFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)

      cs_createNQLFPFiles_channelNum([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, TC)

    
%     %disp(['Creating Spikes day ', num2str(sessionNum)]);
%     %mcz_createNQSpikesFiles([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%     
%     disp(['Creating DIO day ', num2str(sessionNum)]);
%      mcz_createNQDIOFilesFromStateScriptLogs([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)

 end
 
 %% CREATE TETINFO
  createtetinfostruct(dataDir,animID);
 mcz_createcellinfostruct(dataDir,animID,1); %must have spikes
% 
sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodelocation(dataDir,animID,TC,'TC');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');

taskmatrix = [1, 2; 2,2; 3,2; 5,2; 6,2; 7,1; 8,1; 9,1; 9,2; 10,1; 10,2; 10,3; 10,4];
 odorplace = [1, 2; 2,2; 3,2; 5,2; 6,2; 7,1];
 novelodor = [8,1; 9,1];
 novelodor2 = [10,4];
 noodor = [10,1; 10,2; 10,3];
 
cs_riptettag({animID});
riptet = cs_getriptets(animID);

 %% CREATE TASK FILES
%createtaskstruct_nopos(dataDir,'CS41',[1 1; 2 1]);
     createtaskstruct(   dataDir,animID, taskmatrix, 'getcoord_Tmaze'); %index - [day1 epoch1; day1 epoch2; etc]
 
for sessionNum= [1,2,3,5,6,7,8,9,10]%1:numDays
      sj_lineardayprocess(dataDir,animID,sessionNum); %get linpos (need createtaskstruct first)
%      
%       run = taskmatrix((taskmatrix(:,1) == sessionNum), 2)';
%       sj_updatetaskstruct(dataDir,animID,sessionNum, run, 'run'); 
% %      
%       odorplaceepochs = odorplace((odorplace(:,1) == sessionNum), 2)';
%       novelodorepochs = novelodor((novelodor(:,1) == sessionNum), 2)';
%       novelodorepochs2 = novelodor2((novelodor2(:,1) == sessionNum), 2)';
%       noodorepochs = noodor((noodor(:,1) == sessionNum), 2)';
% %      
%      sj_updatetaskenv(   dataDir,animID,sessionNum,odorplaceepochs, 'odorplace');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,novelodorepochs, 'novelodor');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,novelodorepochs2, 'novelodor2');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,noodorepochs, 'noodor');
end

 %% EEG
epochIndex = [];
 for sessionNum= 1:numDays
     cd(topRawDir);
     dayDir = rawDir{sessionNum};
     cd(dayDir)
     recfiles = dir(['*.rec']);
     numEpochs = length(recfiles);
     epochIndex(sessionNum) = numEpochs;
 end
 
refData = {};
for sessionNum = 1:numDays
    dayrefData = repmat(refDataTemplate,epochIndex(sessionNum),1);
    refData{1,sessionNum} = dayrefData;
end

filtPath= ['D:\DataAnalysis\usrlocal\filtering\'];
addpath([dataDir,'EEG']);
for sessionNum= [1,2,3,5:10]%1:numDays
%     disp(['Creating referenced LFP day ', num2str(sessionNum)]);
%     mcz_createRefEEG([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, refData{sessionNum});
%     disp(['Creating theta filtered LFP day ', num2str(sessionNum)]);
%     mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
%     %disp(['Creating delta filtered LFP day ', num2str(sessionNum)]);
%     %mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'deltafilter.mat']);
%     disp(['Creating beta filtered LFP day ', num2str(sessionNum)]);
%     mcz_betadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'betafilter.mat']);
%     disp(['Creating ripple filtered LFP day ', num2str(sessionNum)]);
%     mcz_rippledayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'ripplefilter.mat']);
    disp(['Creating ripple struct day ', num2str(sessionNum)]);
    sj_extractripples(dataDir, animID, sessionNum, riptet, 0.015, 2)
end

%% SPECGRAMS
sjcs_baselinespecgram(animID, 1:numDays, [1,2,3,4], [CA1,PFC], 1,'fpass',[0 40],'movingwin',[1000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, [1,2,3,4], [TC,OB], 1,'fpass',[0 15],'movingwin',[2000 40]/1000);
sjcs_baselinespecgram(animID, 1:numDays, [1,2,3,4], [TC,OB], 1,'fpass',[0 15],'movingwin',[2000 40]/1000,'tapers',[5,1]);

