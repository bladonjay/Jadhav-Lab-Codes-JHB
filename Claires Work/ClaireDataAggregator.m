%% ClaireDataAggregator
%{
this is a description of her dataset and a quick to do list of analyses

the animals:
CS31Expt
CS33Expt
CS34Expt
CS35Expt
CS39Expt
CS41Expt
CS42Expt
CS44Expt

DIOs:
1-	Left reward well
2-	Right reward well
3-	Left Pump
4-	Right Pump
5-	Nosepoke/Buzzer
6-	Left Odor Solenoid
7-	Right Odor Solenoid
8-	Vacuum solenoid
9- sometimes buzzer (itll be an in not an out


The filetypes:
cellinfo: this has general info about the cells:
        spikewidth, mean rate, nspikes, burst prob, tag, area, type
spikes: this corresponds to cellinfo and has:
        data, meanrate, time, x, y, dir, amplitude, posidx, tag, etc.
tetinfo: same, 32 usually, sometimes will have the fields:
        depth, numcells, and area
pos: positional info with fields:
        arg, descript, data (data is time, x, y, dir, vel)
ripples: has a ton of fields: has the ntets organization
        startind,endind,starttime,endtime,midtime,energy,timerange
nosepokewindow: start and stop ts

Rewards:Leftwindows, rightwindows

%}

% first lets take a peek at what our folders have in them:
parentdir='E:\Claire Data';
dirdata=dir(parentdir);
ratfolders=dirdata(contains({dirdata.name},'Expt'));

cellmat={};
for k=1:length(ratfolders)
    cellmat{1,k}=ratfolders(k).name;
    filebits=dir(fullfile(ratfolders(k).folder,ratfolders(k).name));
    filenames={filebits(3:end).name}';
    cellmat(2:(length(filenames)+1),k)=filenames;
end
    % okay this is now a nested struct with all the filenames we'll need

%%

sessnum=1;

for rt=1:length(ratfolders)
    % grab the directory and fill out the name
    thisdir=fullfile(ratfolders(rt).folder, ratfolders(rt).name);
    ratName=ratfolders(rt).name(1:4);
    infofile=[thisdir '\' ratName 'cellinfo.mat'];
    load(infofile); % adds variable 'cellinfo'
    
    % for each session this rat did
    for i=1:length(cellinfo)
        SuperRat(sessnum).name=ratName;
        SuperRat(sessnum).daynum=i;
        SuperRat(sessnum).cellinfo=cellinfo{i};
        SuperRat(sessnum).files.info=infofile;
        fprintf('Now starting rat %s session %d \n',ratName,i);
        
        
        spikefile=[thisdir '\' ratName 'spikes' sprintf('%02d',i) '.mat'];
        if exist(spikefile,'file')==2
            SuperRat(sessnum).files.spikefile=spikefile;
            load(spikefile);
            SuperRat(sessnum).oldspikes=spikes{i};
            % now go in and grab the unit data (across small epochs)...
            cumct=1;
            try
                % for each tetrode in that epoch
                for j=1:length(SuperRat(sessnum).cellinfo{1}) % pull for first epoch first
                    % ok chain goes: epoch, tetrode, unit
                    % i think the same cells are clustered across epochs, so the cell
                    if ~isempty(SuperRat(sessnum).cellinfo{1}{j})
                        % create a set of cells that are accepted
                        for k=1:length(SuperRat(sessnum).cellinfo{1}{j})
                            % the tag means its an accepted cell
                            if isfield(SuperRat(sessnum).cellinfo{1}{j}{k},'tag')
                                % initiate the struct
                                tempstruct=SuperRat(sessnum).cellinfo{1}{j}{k};
                                tempstruct.tet=j;
                                tempstruct.ts=SuperRat(sessnum).oldspikes{1}{j}{k}.data;
                                % now fill in subsequent days
                                for d=2:length(SuperRat(sessnum).oldspikes)
                                    tempstruct.ts=[tempstruct.ts; SuperRat(sessnum).oldspikes{d}{j}{k}.data];
                                end
                                % now a save case cause sometimes its a ref channel
                                if isfield(tempstruct,'descrip')
                                    tempstruct.area=tempstruct.descrip;
                                    tempstruct.type='Other';
                                    tempstruct=rmfield(tempstruct,'descrip');
                                    tempstruct=orderfields(tempstruct,SuperRat(sessnum).units);
                                end
                                SuperRat(sessnum).units(cumct)=tempstruct;
                                cumct=cumct+1;
                            end
                        end
                    end
                end
                fprintf('extracted Units from %s \n',spikefile);
            catch
                fprintf('couldnt extract Units from %s \n',spikefile);
            end
        end
        
        % now go for position data
        posfile=[thisdir '\' ratName 'pos' sprintf('%02d',i) '.mat'];
        if exist(posfile,'file')==2
            SuperRat(sessnum).files.posfile=posfile;
            % load(['E:\Clarie Data\CS31Expt\CS31pos0' num2str(i) '.mat']);
            load(posfile);
            SuperRat(sessnum).oldpos=pos{i};
            % initiate the structure
            try
            temppos=struct('data',SuperRat(sessnum).oldpos{1}.data);
            temppos.data(:,end+1)=1; % this is epoch 1
            temppos.fields={'time','x','y','dir','vel','epoch'};
            temppos.filt(1)=SuperRat(sessnum).oldpos{1}.arg(2);
            temppos.descript=SuperRat(sessnum).oldpos{1}.descript;
            for j=2:length(SuperRat(sessnum).oldpos)
                % gather all the position data
                newpos=SuperRat(sessnum).oldpos{j}.data;
                newpos(:,end+1)=j;
                temppos.data=[temppos.data; newpos];
                temppos.filt(j)=SuperRat(sessnum).oldpos{j}.arg(2);
                temppos.descript(:,j)=SuperRat(sessnum).oldpos{j}.descript;
            end
            SuperRat(sessnum).tracking=temppos;
            SuperRat(sessnum).RunEpochs=2:2:length(SuperRat(sessnum).oldpos);
            fprintf('extracted Coords from %s \n',posfile);
            catch
                fprintf('couldnt extract Coords from %s \n',posfile);
            end
        end
        
        % now for behavior data
        % try a few file formats
        taskfile={[thisdir '\' ratName 'runTrajBounds' sprintf('%02d',i) '.mat']};
        taskfile{2}=[thisdir '\' ratName 'runTrialBounds' sprintf('%02d',i) '.mat'];
        taskfile{3}=[thisdir '\' ratName 'nosepokeWindow' sprintf('%02d',i) '.mat'];
        taskfile{4}=[thisdir '\' ratName 'rewardTimes' sprintf('%02d',i) '.mat'];
        taskfile{5}=[thisdir '\' ratName 'rewards' sprintf('%02d',i) '.mat'];
        taskfile{5}=[thisdir '\' ratName 'odorTriggers' sprintf('%02d',i) '.mat'];
        
        % put this into its own fx because its a lot of accounting
        [trialdata,usedfiles,oldbeh] = ParseClaireBehavior(taskfile,i);
        if ~isempty(usedfiles)
            SuperRat(sessnum).files.taskfile=usedfiles;
            SuperRat(sessnum).trialdata=trialdata;
            SuperRat(sessnum).oldbeh=oldbeh;
            if iscell(usedfiles)
                for k=1:length(usedfiles)
                    fprintf('extracted behavior data from %s \n', usedfiles{k});
                end
            else
                fprintf('extracted behavior data from %s \n', usedfiles);
            end
        else
            fprintf('None of the behavior files could be found \n');
        end
        
        % and finally get the ripples
        ripfile=[thisdir '\' ratName 'ripples' sprintf('%02d',i) '.mat'];
        if exist(ripfile,'file')==2
            SuperRat(sessnum).files.ripfile=ripfile;
            load(ripfile);
            % now go and concatenate the ripple times for all the tetrodes
            % grab that day
            % need to collapse each session into a single tetrode
            tetcat=1;
            for tet=1:length(ripples{i}{1})
                if ~isempty(ripples{i}{1}{tet})
                    ripples{i}{1}{tet}.tetrode=tet;
                    if tetcat==1, ripdata=ripples{i}{1}{tet}; else, ripdata(tetcat)=ripples{i}{1}{tet}; end 
                    
                    for ses=2:length(ripples{i})
                        % i actually dont even know what half of these
                        % are...
                        ripdata(tetcat).starttime=[ripdata(tetcat).starttime; ripples{i}{ses}{tet}.starttime];
                        ripdata(tetcat).endtime=[ripdata(tetcat).endtime; ripples{i}{ses}{tet}.endtime];
                        ripdata(tetcat).midtime=[ripdata(tetcat).midtime; ripples{i}{ses}{tet}.midtime];
                        ripdata(tetcat).peak=[ripdata(tetcat).peak; ripples{i}{ses}{tet}.peak];
                        ripdata(tetcat).energy=[ripdata(tetcat).energy; ripples{i}{ses}{tet}.energy];
                        ripdata(tetcat).maxthresh=[ripdata(tetcat).maxthresh; ripples{i}{ses}{tet}.maxthresh];
                    end
                    tetcat=tetcat+1;
                end
            end
            SuperRat(sessnum).ripdata=ripdata;
            fprintf('extracted ripple data from %s \n', ripfile);
        end
        % and of course advance the session number
        sessnum=sessnum+1;
        fprintf('finished Rat %s %d \n \n',ratName,i);
    end
end
% now go in and grab the unit data (across small epochs)...

    
    
%% this is because i had to add in behavior data for CS44 at the very end


sessnum=21; % starting session number
ratName=SuperRat(sessnum).name;
thisdir=SuperRat(sessnum).files.info;
thisdir=thisdir(1:find(thisdir=='\',1,'last')-1);

infofile=[thisdir '\' ratName 'cellinfo.mat'];
load(infofile); % adds variable 'cellinfo'
oldbeh={};
% for each session this rat did
for i=1:length(cellinfo)
    
    taskfile={[thisdir '\' ratName 'runTrajBounds' sprintf('%02d',i) '.mat']};
    taskfile{2}=[thisdir '\' ratName 'runTrialBounds' sprintf('%02d',i) '.mat'];
    taskfile{3}=[thisdir '\' ratName 'nosepokeWindow' sprintf('%02d',i) '.mat'];
    taskfile{4}=[thisdir '\' ratName 'rewardTimes' sprintf('%02d',i) '.mat'];
    taskfile{5}=[thisdir '\' ratName 'rewards' sprintf('%02d',i) '.mat'];
    taskfile{6}=[thisdir '\' ratName 'odorTriggers' sprintf('%02d',i) '.mat'];
    
    
    % put this into its own fx because its a lot of accounting
    [trialdata,usedfiles,thisoldbeh] = ParseClaireBehavior(taskfile,i);
    
    if ~isempty(usedfiles)
        SuperRat(sessnum).files.taskfile=usedfiles;
        SuperRat(sessnum).trialdata=trialdata;
        SuperRat(sessnum).oldbeh=thisoldbeh(:,end);
        if iscell(usedfiles)
            for k=1:length(usedfiles)
                fprintf('extracted behavior data from %s \n', usedfiles{k});
            end
        else
            fprintf('extracted behavior data from %s \n', usedfiles);
        end
    else
        fprintf('None of the behavior files could be found \n');
    end
    sessnum=sessnum+1;
end


%% this is because CS39 doesnt have the reward10 flag, so I have to make my own\
%{
basically you can tell if there is a reward based on whether the
immediately following event is a reward, and the reward says which side
it is.
so there are 96 sniffstart and sniffends
there are 69 rewardstarts and rewardends
the left/right apply to the rewards, and obviously if they're rewards,


The rewards are only for correct trials, so only when the reward was actually
dispensed and not when the animal poked at the reward but didn't get anything.
So in that case, the r/l rewards correspond to r/l odors, and in trials
that dont follow with a reward, the RL odor ID is unknown** (which is okay
because you really only want to use correct trials anyways
%}

for ses=1:length(SuperRat)
    % if we dont have good data
    if ~isfield(SuperRat(ses).trialdata,'CorrIncorr10')
        rawtrial=SuperRat(ses).trialdata;
        SuperRat(ses).trialdata.starttime=rawtrial.sniffstart;
        SuperRat(ses).trialdata.endtime=rawtrial.sniffend;
        
        % now have to extract left and right from old data
        oldbeh=SuperRat(ses).oldbeh(:,end);
        lefts=[]; rights=[];
        allbeh=[oldbeh{:}];
        for k=1:length(allbeh)
            if isstruct(allbeh{k})
            lefts=[lefts; allbeh{k}.leftWindows];
            rights=[rights; allbeh{k}.rightWindows];
            end
        end
        lefts=unique(lefts,'rows'); rights=unique(rights,'rows');
        allcorrect=sortrows([lefts ones(length(lefts),1); rights zeros(length(rights),1)]);
        % now go through each correct and find the most recent study and
        % asign it a 1
        leftright10=nan(length(SuperRat(ses).trialdata.sniffstart),1);
        CorrIncorr10=zeros(length(SuperRat(ses).trialdata.sniffstart),1);

        for k=1:length(allcorrect)
            lastodor=find(SuperRat(ses).trialdata.sniffstart<allcorrect(k,1),1,'last');
            leftright10(lastodor)=allcorrect(k,3);
            CorrIncorr10(lastodor)=1;
        end
        SuperRat(ses).trialdata.leftright10=leftright10;
        SuperRat(ses).trialdata.CorrIncorr10=CorrIncorr10;
        fprintf('fixed beh from %s Ses %d, %d trials %d correct \n',...
            SuperRat(ses).name,ses,length(SuperRat(ses).trialdata.CorrIncorr10),...
        sum(SuperRat(ses).trialdata.CorrIncorr10));
    end        
end

%% need to designate the epoch for all the sessions

for i=1:length(SuperRat)
    tsdata=SuperRat(i).tracking.data;
    EpochChangeInds=find(diff(tsdata(:,6))~=0);
    Trialts=SuperRat(i).trialdata.sniffstart;
    Trialts(:,2)=1;
    for g=1:length(EpochChangeInds)
        Trialts(Trialts(:,1)>=tsdata(EpochChangeInds(g)+1,1),2)=g+1;
    end
    Trialts(Trialts(:,1)>=tsdata(end,1),2)=nan;
    SuperRat(i).trialdata.EpochInds=Trialts;
end





















