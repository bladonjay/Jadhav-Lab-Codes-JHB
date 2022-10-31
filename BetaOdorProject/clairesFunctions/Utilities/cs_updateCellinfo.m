%cs_updateCellinfo

topDir = cs_setPaths();

%% CS35
animID = 'CS35';

topDir = cs_setPaths();
topRawDir = [topDir,animID,'Expt\',animID,'\'];
dataDir = [topDir,animID,'Expt\',animID,'_direct\'];

CA1   = [1 2 3 4 17 18 20 27 28 29 30 31 32];
 hpcRef=  19;
 PFC   = [7 8 10 11 13 14 15 16 21 23 24 25 26];
 pfcRef=   22;
 OB = [12];
 obRef = 12;
 
 
mcz_createcellinfostruct(dataDir,animID); %must have spikes
createtetinfostruct(dataDir,animID);

sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');

cs_addNumCells(dataDir,animID)

%% CS39
animID= 'CS39';     % your animals name

topDir = cs_setPaths();
topRawDir = [topDir,animID,'Expt\',animID,'\'];
dataDir = [topDir,animID,'Expt\',animID,'_direct\'];

% tetrode metadata (which tetrode is in which area)
CA1   = [16 17 19 20 27 28 29 30 4 5];
 hpcRef=  31;
 PFC   = [7 8 9 10 21 26 11 12 14 15];
 pfcRef=   13;
 OB = [22 23 24 25];
 obRef = 22;
 V1 = [1 2 3 32];
 v1Ref = 1;
 
mcz_createcellinfostruct(dataDir,animID); %must have spikes
createtetinfostruct(dataDir,animID);

sj_addtetrodelocation(dataDir,animID,OB,'V1');
sj_addtetrodelocation(dataDir,animID,OB,'v1Ref');
sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');

cs_riptettag({animID});
 
%% CS41
animID= 'CS41';    % your animal's name
topRawDir= [topDir,animID,'Expt\',animID,'\']; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir = [topDir,animID,'Expt\',animID, '_direct\'];

% tetrode metadata (which tetrode is in which area)
TC = 24; 
OB = 12;
CA1 = [17 18 19 20 28 29 30 32 1 2 3 4];
hpcRef = 27;
PFC = [5 7 8 10 21 22 23 25 26 11 13 14 15 16];
pfcRef = 9;

mcz_createcellinfostruct(dataDir,animID,1); %must have spikes
%createtetinfostruct(dataDir,animID);

sj_addtetrodelocation(dataDir,animID,OB,'TC');
sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');


cs_riptettag({animID});


%% CS42
animID= 'CS42';    % your animal's name
topRawDir= [topDir,animID,'Expt\',animID,'\']; % raw directory (contains one folder for each day of an experiment ONLY)
dataDir = [topDir,animID,'Expt\',animID, '_direct\'];

% tetrode metadata (which tetrode is in which area)
TC = 24; 
OB = 12;
CA1 = [17 18 20 27 28 29 30 32 2 3 4];
hpcRef = 1;
PFC = [6 7 9 10 21 22 23 25 26 11 13 14 16];
pfcRef = 8;

mcz_createcellinfostruct(dataDir,animID,1); 
createtetinfostruct(dataDir,animID);

sj_addtetrodelocation(dataDir,animID,OB,'TC');
sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');

cs_riptettag({animID});

%% CS44

animID = 'CS44';
topRawDir= [topDir,animID,'Expt\',animID,'\']; 
dataDir = [topDir,animID,'Expt\',animID, '_direct\'];

OB = 11;
CA1 = [17 18 19 20 27 28 30 31 32 1 2 3 4];
hpcRef = 29;
PFC = [5 7 8 9 10 21 22 23 25 26 12 13 14 15];
pfcRef = 24;

mcz_createcellinfostruct(dataDir,animID,1); %must have spikes
createtetinfostruct(dataDir,animID);

sj_addtetrodelocation(dataDir,animID,OB,'OB');
sj_addtetrodelocation(dataDir,animID,CA1,'CA1');
sj_addtetrodelocation(dataDir,animID,PFC,'PFC');
sj_addtetrodedescription(dataDir,animID,hpcRef,'hpcRef');
sj_addtetrodedescription(dataDir,animID,pfcRef,'pfcRef');

cs_riptettag({animID});

%% All animals
cs_findBetaTets({'CS35'},'odorplace');
cs_listAllCells
cs_cellTypeTag
cs_listNPCells
cs_cellSelectivityTag
cs_listSelectiveCells