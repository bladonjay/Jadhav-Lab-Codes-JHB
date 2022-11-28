%% METADATA

% have Trodes, usrlocal, Src_matlab, TrodesToMatlab in path
animID= 'CS44';    % your animal's name

topDir = cs_setPaths(); %set this manually if not on Claire's computer

filtPath= 'D:\DataAnalysis\usrlocal\filtering\'; %set manually - path to usrlocal\filtering folder

%topRawDir= [topDir,animID,'Expt\',animID,'\']; % raw directory (contains one folder for each day of an experiment ONLY)
topRawDir = 'D:\OdorPlaceAssociation\RAW\CS44\';
dataDir = [topDir,animID,'Expt\',animID, '_direct\'];

% tetrode metadata (which tetrode is in which area) 
OB = 11;
CA1 = [17 18 19 20 27 28 30 31 32 1 2 3 4];
hpcRef = 29;
PFC = [5 7 8 9 10 21 22 23 25 26 12 13 14 15];
pfcRef = 24;

%Create refDataTemplate: array indicating which tetrode should be used as
%reference. 0 indicates tetrode was missing, or was a reference.
allTets = [OB,CA1,hpcRef,PFC,pfcRef];
refDataTemplate = 1:32;
unused = find(~ismember(refDataTemplate,allTets));
refDataTemplate([OB, hpcRef, pfcRef,unused]) = 0;
refDataTemplate(CA1) = hpcRef;
refDataTemplate(PFC) = pfcRef;


%taskmatrix = [1, 2; 1,4; 2,2; 2,4; 3,2; 3,4; 4,2; 5,2;  6,2; 6,4; 6,6; 7,2; 7,4; 7,6; 8,2; 8,4; 9,2; 9,4; 9,6; 10,2; 10,4; 10,6];
taskmatrix = [1, 2; 1,4; 2,2; 2,4; 3,2; 3,4; 4,2; 5,2; 6,2; 6,4; 6,6; 7,2; 7,4];
odorplace = [1, 2; 1,4; 2,2; 2,4; 3,2; 3,4; 4,2; 5,2; 8,2; 10,6];
 reversefam = [10,2; 10,4];
 novelodor = [6,2; 6,4; 6,6];
 novelodor2 = [7,2; 7,4];
 noodor = [8,4];
 reversenovel1 = [9,6];
 reversenovel2 = [9,2; 9,4];
 novelmix = [7,6];
 


 % experiment metadata
 cd(topRawDir);
 rawDir=dir();
    rawDir= {rawDir(3:end).name};
 numDays= length(rawDir);
 % The above code assumes raw directory contains one folder for each day of experiment ONLY
 
 
 %% DAY DEPENDENT - create eeg files, and filtered lfp files
  for sessionNum= 8%1:numDays
    
%     disp(['Creating pos day ', num2str(sessionNum)]); 
%     mcz_createNQPosFiles_new([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, [])
%     
    disp(['Creating LFP day ', num2str(sessionNum)]);
    mcz_createNQLFPFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%     
% %     disp(['Creating Spikes day ', num2str(sessionNum)]);
% %     mcz_createNQSpikesFiles([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
% %     
    disp(['Creating DIO day ', num2str(sessionNum)]);
    mcz_createNQDIOFilesFromStateScriptLogs([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)

  end
 
%  
 %createtaskstruct(   dataDir,animID, taskmatrix, 'getcoord_Tmaze'); %index - [day1 epoch1; day1 epoch2; etc]
 
 
 %% CREATE TETINFO
 createtetinfostruct(dataDir,animID,1);
%mcz_createcellinfostruct(dataDir,animID); %must have spikes

sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');



 %% CREATE TASK FILES
createtaskstruct_nopos(dataDir,'CS44',noodor);
     
for sessionNum= 8
      %sj_lineardayprocess(dataDir,animID,sessionNum) %get linpos (need createtaskstruct first)
%      
%       run = taskmatrix((taskmatrix(:,1) == sessionNum), 2)';
%       sj_updatetaskstruct(dataDir,animID,sessionNum, run, 'run'); 
% %      
%       odorplaceepochs = odorplace((odorplace(:,1) == sessionNum), 2)';
%       novelodorepochs = novelodor((novelodor(:,1) == sessionNum), 2)';
%       novelodorepochs2 = novelodor2((novelodor2(:,1) == sessionNum), 2)';
      noodorepochs = noodor((noodor(:,1) == sessionNum), 2)';
%       reversefamepochs = reversefam((reversefam(:,1) == sessionNum), 2)';
%       reversenov1epochs = reversenovel1((reversenovel1(:,1) == sessionNum), 2)';
%       reversenov2epochs = reversenovel2((reversenovel2(:,1) == sessionNum), 2)';
%       novelmixepochs = novelmix((novelmix(:,1) == sessionNum), 2)';
%      
%      sj_updatetaskenv(   dataDir,animID,sessionNum,odorplaceepochs, 'odorplace');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,novelodorepochs, 'novelodor');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,novelodorepochs2, 'novelodor2');
     sj_updatetaskenv(   dataDir,animID,sessionNum,noodorepochs, 'noodor');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,noodorepochs, 'reversefam');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,noodorepochs, 'reversenovel1');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,noodorepochs, 'reversenovel2');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,noodorepochs, 'novelmix');
     
     %cs_addPosToSpikesFiles(animID, dataDir, sessionNum)
end

 %% EEG
epochIndex = [];
 for sessionNum= 8%1:numDays
     cd(topRawDir);
     dayDir = rawDir{sessionNum};
     cd(dayDir)
     recfiles = dir(['*.rec']);
     numEpochs = length(recfiles);
     epochIndex(sessionNum) = numEpochs;
 end

 
refData = {};
for sessionNum = 8%1:numDays
    dayrefData = repmat(refDataTemplate,epochIndex(sessionNum),1);
    refData{1,sessionNum} = dayrefData;
end

filtPath= ['D:\DataAnalysis\usrlocal\filtering\'];
addpath([dataDir,'EEG']);
riptet = cs_getriptets(animID);

 for sessionNum= 8%1:numDays
    disp(['Creating filtered LFP day ', num2str(sessionNum)]);
    mcz_createRefEEG([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, refData{sessionNum});
%     disp(['Creating theta filtered LFP day ', num2str(sessionNum)]);
%     mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
%     %mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'deltafilter.mat']);
    disp(['Creating beta filtered LFP day ', num2str(sessionNum)]);
     mcz_betadayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'betafilter.mat']);
    disp(['Creating ripple filtered LFP day ', num2str(sessionNum)]);
    mcz_rippledayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'ripplefilter.mat']);
     %sj_extractripples(dataDir, animID, sessionNum, riptet, 0.015, 2)
 end

%% SPECGRAMS
sjcs_baselinespecgram(animID, 1:5, [2,4], [CA1,PFC,OB], 1,'fpass',[0 40],'movingwin',[1000 20]/1000);
%sjcs_baselinespecgram(animID, 1:5, [2,4], [OB], 1,'fpass',[0 15],'movingwin',[2000 40]/1000);
