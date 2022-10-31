%% METADATA

% have Trodes, usrlocal, Src_matlab, matclust, TrodesToMatlab, and Pipeline in path
% To extract matclust files for clustering, run createAllMatclustFiles.m (located in TrodesToMatlab)

% path and animal metadata
topRawDir= 'D:\OdorPlaceAssociation\CS35Expt\CS35\'; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir= 'D:\OdorPlaceAssociation\CS35Expt\CS35_direct\';    % data directory
animID= 'CS35';     % your animals name

% tetrode metadata (which tetrode is in which area, which is reference)
 CA1   = [1 2 3 4 17 18 20 27 28 29 30 31 32];
 hpcRef=  19;
 PFC   = [7 8 10 11 13 14 15 16 21 23 24 25 26];
 pfcRef=   22;
 OB = [12];
 obRef = 12;
 
 % experiment metadata
 cd(topRawDir);
 rawDir=dir();
    rawDir= {rawDir(3:end).name};
 numDays= length(rawDir);
 % The above code assumes raw directory contains one folder for each day of experiment ONLY

 
     
     
 
 %% DAY DEPENDENT NON EEG- pos, lfp, spikes, dio, task, linpos
 for sessionNum= 1:numDays;
     
     
     % CONVERT RAW DATA TO BASIC FF DATA
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cs_addPosToSpikesFiles(animID, dataDir, sessionNum)
     %mcz_createNQPosFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, [], 'diodenum', diodenum)

    % TRANSITION TO THIS WHEN USING RANGE DURING POSITION TRACKING INSTEAD OF VIDEO+FIRST FRAME METHOD
     %mcz_createNQPosFiles_new([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, [])
% 
      %mcz_createNQLFPFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%     mcz_createNQSpikesFiles([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
     %mcz_createNQDIOFiles(    [topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
%      mcz_createNQDIOFilesFromStateScriptLogs([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum)
     
     
     
     % FF SECONDARY DATA STRUCTURES
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %riptet = [18,20,31];
     % sj_extractripples(dataDir, animID, sessionNum, riptet, 0.015, 2)

 end
 
%  mcz_createNQSpikesFiles([topRawDir rawDir{3}], dataDir, animID, sessionNum)
 % FF BEHAVIOR METADATA STRUCTURES
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      createtaskstruct(   dataDir,animID, [1 2; 1 4; 2 2; 2 4; 3 2; 3 4], 'getcoord_Tmaze'); %index - [day1 epoch1; day1 epoch2; etc]
%     createtaskstruct_nopos( dataDir,animID, [1 2; 1 4; 2 2; 2 4; 3 2; 3 4]); 
for sessionNum= 1:3 
     % sj_lineardayprocess(dataDir,animID,sessionNum) %get linpos (need createtaskstruct first)
     %sj_updatetaskstruct(dataDir,animID,sessionNum,[2,4], 'run'); % one
%      %line for each day - need to fix so it works regardless of
%      % number of epochs per day
%      sj_updatetaskenv(   dataDir,animID,sessionNum,[2,4], 'odorplace');
 

%sj_extractripples(dataDir, animID, sessionNum, riptet, 0.015, 2)
end


     
     
     %sj_updatetaskenv(   dataDir,animID,sessionNum,[#],
     %'wtr2');
 %% DAY INDEPENDENT- cellinfo, tetinfo
 
%  odorPlaceTaskDays = [1,2,3];
% save([dataDir,animID,'odorPlaceTaskDays.mat'],'odorPlaceTaskDays'); %change numbers in matrix to correct days

 
 % FF CELL AND TET METADATA STRUCTURES
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  createtetinfostruct(dataDir,animID);
  mcz_createcellinfostruct(dataDir,animID); %must have spikes
%  cs_cellTypeTag
%  cs_cellSelectivityTag
 cs_addNumCells(dataDir,animID);
%  % CELL AND TET METADATA POPULATION
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %riptet = [18,20, 28,32];
 
%sj_addtetrodedescription(dataDir,animID,riptet,'riptet');
%cs_addriptet(dataDir,animID)
      sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
     sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');
       sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
      sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
       sj_addtetrodelocation(dataDir,animID,OB,'OB');
      sj_addtetrodedescription(dataDir,animID,obRef,'obRef');
%       betatets = [8, 11, 12, 14, 16, 18, 20, 24, 26, 28, 31, 32];
% cssj_addtetrodedescription2(dataDir,'CS35',betatets,'betatet');
% %   
%  sj_addtetrodedescription(dataDir,animID,Ref,'Ref');
 %sj_addtetrodedescription(dataDir,animID,pfcRef,'PFCRef');
 %sj_addcellinfotag2(dataDir,animID);
 
%sj_extractripples(dataDir, animID, sessionNum, riptet, 0.015, 2)

 %% EEG PRIMARY PIPELINE
 
%refData -- an EX32 matrix with the local reference for each tetrode, where
%unused or reference tetrodes have a ref of zero. E is 1:total # of epochs
% TODO- add changes over days? 
 
% filtering/filter path

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
 
refDataTemplate = [19 19 19 19 0 0 22 22 22 22 22 0 22 22 22 22 19 19 0 19 22 0 22 22 22 22 19 19 19 19 19 19];
refData = {};
for sessionNum = 1:numDays
    dayrefData = repmat(refDataTemplate,epochIndex(sessionNum),1);
    refData{1,sessionNum} = dayrefData;
end

filtPath= ['D:\DataAnalysis\usrlocal\filtering\'];
addpath('D:\OdorPlaceAssociation\CS35Expt\CS35_direct\EEG');
for sessionNum= 1:numDays

% mcz_createRefEEG([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, refData{sessionNum});
%      % EEG wrt ground only
%mcz_deltadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'deltafilter.mat']);
%  disp(['doing day ',num2str(sessionNum), ' theta']);
%  mcz_thetadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'thetafilter.mat']);
% 
% disp(['doing day ',num2str(sessionNum), ' beta']);
% mcz_betadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'betafilter.mat']);

%      % EEG wrt ground and local reference
% disp(['doing day ',num2str(sessionNum), ' gamma']);
% mcz_gammadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'gammafilter.mat']);
% % 
% disp(['doing day ',num2str(sessionNum), ' ripple']);
% mcz_rippledayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'ripplefilter.mat']);

% disp(['doing day ',num2str(sessionNum), ' betaRef']);
% mcz_betadayprocess_Ref([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'betafilter.mat']);


%mczcs_lowgammadayprocess([topRawDir rawDir{sessionNum}], dataDir, animID, sessionNum, 'f', [filtPath 'lowgammafilter.mat']);
end
% 
sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB], 1,'fpass',[0 10],'movingwin',[2000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB], 1,'fpass',[0 40],'movingwin',[1000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB, hpcRef, pfcRef], 0,'fpass',[0 40],'movingwin',[1000 20]/1000);
sjcs_baselinespecgram(animID, 1:numDays, [2 4], [CA1, PFC, OB, hpcRef, pfcRef], 0,'fpass',[0 100],'movingwin',[400 20]/1000);

