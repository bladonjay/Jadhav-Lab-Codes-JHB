%% METADATA

% have Trodes, usrlocal, Src_matlab, matclust, TrodesToMatlab, and Pipeline in path
% To extract matclust files for clustering, run createAllMatclustFiles.m (located in TrodesToMatlab)

% path and animal metadata
animID= 'CS39';     % your animals name

topDir = cs_setPaths();
topRawDir = [topDir,animID,'Expt\',animID,'\'];
dataDir = [topDir,animID,'Expt\',animID,'_direct\'];

% tetrode metadata (which tetrode is in which area, which is reference)
 CA1   = [16 17 19 20 27 28 29 30 4 5];
 hpcRef=  31;
 PFC   = [7 8 9 10 21 26 11 12 14 15];
 pfcRef=   13;
 OB = [22 23 24 25];
 obRef = 22;
 V1 = [1 2 3 32];
 v1Ref = 1;
 
 %riptet = [16, 17, 19, 20, 27, 29, 30, 4, 5];
 
 refDataTemplate = [0 1 1 31 31 0 13 13 13 13 13 13 0 13 13 31 31 0 31 31 13 0 22 22 22 13 31 31 31 31 0 1];

 
 % experiment metadata
 cd(topRawDir);
 rawDir=dir();
    rawDir= {rawDir(3:end).name};
 numDays= length(rawDir);
 % The above code assumes raw directory contains one folder for each day of experiment ONLY

 taskmatrix = [1, 2; 1,4; 2,2; 2,3;  3,2; 3,3; 3,5; 4,2; 4,3; 4,5; 5,2; 5,3; 5,5; 6,2; 6,3; 6,5; 7,2; 7,4; 7,5; 7,6; 7,8];
 odorplace = [1,2; 1,4; 2,2; 3,2; 3,5; 4,2; 5,2; 6,2; 7,2 ; 7,4];
 novelodor = [2,3; 3,3; 4,3; 4,5; 5,3; 5,5; 6,3; 6,5; 7,5; 7,6];
 noodor = [7,8];
 
 %% DAY DEPENDENT NON EEG- pos, lfp, spikes, dio, task, linpos
 for sessionNum= 1:numDays
    cs_addPosToSpikesFiles(animID, dataDir, sessionNum)
%     disp(['Creating pos day ', num2str(sessionNum)]); 
%       mcz_createNQPosFiles_new([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, [])
%     
%     disp(['Creating LFP day ', num2str(sessionNum)]);
%     mcz_createNQLFPFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%     
%     disp(['Creating Spikes day ', num2str(sessionNum)]);
%     %mcz_createNQSpikesFiles([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%     
%     disp(['Creating DIO day ', num2str(sessionNum)]);
    % mcz_createNQDIOFilesFromStateScriptLogs([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%      
%      
%      
  
 end
 
     % FF BEHAVIOR METADATA STRUCTURES
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     createtaskstruct(   dataDir,animID, [7,2; 7,4; 7,5; 7,6; 7,8], 'getcoord_Tmaze'); %index - [day1 epoch1; day1 epoch2; etc]
 for sessionNum= 1:numDays
      %sj_lineardayprocess(dataDir,animID,sessionNum) %get linpos (need createtaskstruct first)
%      
%      run = taskmatrix((taskmatrix(:,1) == sessionNum), 2)';
%      sj_updatetaskstruct(dataDir,animID,sessionNum, run, 'run'); 
%      
%      odorplaceepochs = odorplace((odorplace(:,1) == sessionNum), 2)';
%      novelodorepochs = novelodor((novelodor(:,1) == sessionNum), 2)';
      noodorepochs = noodor((noodor(:,1) == sessionNum), 2)';
%      
%      sj_updatetaskenv(   dataDir,animID,sessionNum,odorplaceepochs, 'odorplace');
%      sj_updatetaskenv(   dataDir,animID,sessionNum,novelodorepochs, 'novelodor');
      sj_updatetaskenv(   dataDir,animID,sessionNum,noodorepochs, 'noodor');
 end
%      
%  %% DAY INDEPENDENT- cellinfo, tetinfo
%  
%  % FF CELL AND TET METADATA STRUCTURES
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  createtetinfostruct(dataDir,animID);
  mcz_createcellinfostruct(dataDir,animID); %must have spikes
 cs_addNumCells(dataDir,animID);
%  
%  % CELL AND TET METADATA POPULATION
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   riptet = [16, 17, 19, 20, 27, 29, 30, 4, 5];
% %  
    %sj_addtetrodedescription(dataDir,animID,riptet,'riptet');
    sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
    sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');
    sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
    sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
    sj_addtetrodelocation(dataDir,animID,OB,'OB');
    sj_addtetrodedescription(dataDir,animID,obRef,'obRef');
      sj_addtetrodelocation(dataDir,animID,V1,'V1');
     sj_addtetrodedescription(dataDir,animID,v1Ref,'v1Ref');
% 
cs_riptettag({animID});
riptets = cs_getriptets(animID);
% %  betatets = [8 10 18 20 22 24 25 30 31];
% %  cssj_addtetrodedescription2(dataDir,'CS33',betatets,'betatet');


 %% EEG PRIMARY PIPELINE
 
%Creates refData array with different numbers of copies of refDataTemplate
%on different days according to number of epochs- siplifies so script
%doesn't have to be run multiple times for different days. 


%Create epoch index (number of epochs per day)
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
addpath('D:\OdorPlaceAssociation\CS39Expt\CS39_direct\EEG');
for sessionNum= 1:numDays
disp(['doing day ',num2str(sessionNum),]);
    
% mcz_createRefEEG([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, refData{sessionNum});
% 
% % EEG wrt ground only
% mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'deltafilter.mat']);
% mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
% % disp(['doing day ',num2str(sessionNum), ' beta']);
% mcz_betadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'betafilter.mat']);
% 
% % EEG wrt ground and local reference
% mcz_gammadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'gammafilter.mat']);
% mcz_rippledayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'ripplefilter.mat']);
     sj_extractripples(dataDir, animID, sessionNum, riptets, 0.015, 2)

end
% 
sjcs_baselinespecgram(animID, 1:numDays, 2, [CA1, PFC, OB, V1], 1,'fpass',[0 20],'movingwin',[2000 40]/1000);
sjcs_baselinespecgram(animID, 1:numDays, 2, [CA1, PFC, OB, V1], 1,'fpass',[0 40],'movingwin',[1000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, 2, [CA1, PFC, OB, V1, hpcRef, pfcRef, obRef, v1Ref], 0,'fpass',[0 40],'movingwin',[1000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, 2, [CA1, PFC, OB, V1, hpcRef, pfcRef, obRef, v1Ref], 0,'fpass',[0 100],'movingwin',[400 20]/1000);

sjcs_baselinespecgram(animID, 1:numDays, 2, [V1], 1,'fpass',[0 20],'movingwin',[2000 40]/1000);
sjcs_baselinespecgram(animID, 1:numDays, 2, [V1], 1,'fpass',[0 40],'movingwin',[1000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, 2, [obRef, v1Ref], 0,'fpass',[0 40],'movingwin',[1000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, 2, [obRef, v1Ref], 0,'fpass',[0 100],'movingwin',[400 20]/1000);

