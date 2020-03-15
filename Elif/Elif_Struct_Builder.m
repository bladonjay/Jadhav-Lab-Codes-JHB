% Elif data workup

% First load up the experimental metadata
load('E:\Elif DG GABAa project\Animal Data\animal_metadata_190421.mat')


% hardcoded folders
parentdir='E:\Elif DG GABAa project\Animal Data';

% go into it and grab the files
folders=dir(parentdir);%
%
mymice=folders(contains({folders.name},'Exp'));
% now cycle through and grab the data

mousefiles=struct('directory',[]);
for k=1:length(mymice)
    % first get the path
    mousefiles(k).directory=fullfile(mymice(k).folder,mymice(k).name);
    % pull all filenames
    filebits=getAllFiles(fullfile(mymice(k).folder,mymice(k).name),1);
    % now that we have the directory, just get the file paths in that
    % directory
    filebits=cellfun(@(a) a(length(mousefiles(k).directory)+2:end), filebits, 'UniformOutput', false);
    % now sort them
    %  tracking
    mousefiles(k).tracking=filebits(cellfun(@(a) contains(a,'pos'), filebits));
    % and finally, spikes
    mousefiles(k).unitdata=filebits(cellfun(@(a) contains(a,'spike'), filebits));
    % behavioral data
    mousefiles(k).taskdata=filebits(cellfun(@(a) contains(a,'task'), filebits));
    % and all eegdata
    mousefiles(k).eeg=filebits(cellfun(@(a) contains(a,'eeg'), filebits));
    % ripple data
    mousefiles(k).rippledata=filebits(cellfun(@(a) contains(a,'rip'), filebits));
end

%%  this is the metadata for most of the animals

% this is a struct with the following organization:
%{
animDB                  -struct
    animal              -string
    genotype            -string
    gender              -string
    animal uid          -string
    DOB                 -datetime
    Project             -string
    experiment_dir      -string
    preeimplant_weight  -string or double
    implant_date        -datetime
    implant_age         -double         (days?)
    implant_weight      -string or double

    recording_data      -struct
        day             -double
        date            -datetime
        time            -datetime
        weight          -double or empty
        age             -double         (days?)
        estrus          -double

        epochs          -struct
            epoch       -double
            epoch_type  -string
            environment -string
            comments    -string or empty

        tet_info        -struct
            tetrode
            depth
            target
            riptet      -logical
            ref         -logical
            single_unit -logical
            multi_unit  -logical
            exclude     -logical


so basically there is more information in this struct than you would ever
need to know, and you can probably use this as a reference to verify all
the following datafiles in each mouse folder

%}
%%
% it looks like the animals have great metadata, so i'll basically reform
% that struct and add in spike and lfp data when i want to

sesstally=1;
clear SuperMouse;
for ms=1:length(animDB)
    % create a struct where we drop metadata into one field
    thismouse=struct('name',animDB(ms).animal,'mouse_meta',rmfield(animDB(ms),{'animal','recording_data'}));
    % now we can go through each recording
    for day=1:length(animDB(ms).recording_data)
        % reorganize the fields in this struct
        thisDay=thismouse;
        TDBindex=find(cellfun(@(a) contains(a,'EE8'), {mousefiles.directory}));
        thisDay.day_meta=rmfield(animDB(ms).recording_data(day),{'epochs','tet_info','date'});
        
        daynum=sprintf('%02d',thisDay.day_meta.day);
        thisDay.daynum=daynum;
        thisDay.day=string(animDB(ms).recording_data(day).date);
        thisDay.epochs=animDB(ms).recording_data(day).epochs;
        thisDay.tetinfo=animDB(ms).recording_data(day).tet_info;
        
        % add the task types that occurred today
        thisDay.taskTypes=unique({thisDay.epochs.epoch_type});
        
        % and prefill pos tracking and unit data
        thisDay.unitdata=[];
        thisDay.oldspikes=[];
        thisDay.rawtracking=[];
        thisDay.coords=[];
        thisDay.sourcefiles=[];
        thisDay.lincoords=[];
        thisDay.rawlintracking=[];
        % now we have the chance to add tracking and unit data.
        % the unit data will be for the day, whereas the
        % tracking data will follow the epochs
        % we'll use the day (turn into 2 digit with leading zero)
        
        % search for the pos data (without lin)
        msindex=find(cellfun(@(a) contains(a,thismouse.name), {mousefiles.directory}));
        posfile=find(contains(mousefiles(msindex).tracking,['pos' daynum]) & ...
            ~contains(mousefiles(msindex).tracking,'lin'));
        if ~isempty(posfile)
            rawtracking=load(fullfile(mousefiles(msindex).directory,mousefiles(msindex).tracking{posfile}));
            thisDay.rawtracking=rawtracking.pos(~isempty(rawtracking.pos));
            thisDay.sourcefiles=[thisDay.sourcefiles; {fullfile(mousefiles(msindex).directory,mousefiles(msindex).tracking{posfile})}];
            epochpos=rawtracking.pos{cellfun(@(a) ~isempty(a),rawtracking.pos)};
            
            clear rawtracking;
            for k=1:length(epochpos)
                rawtracking(k)=epochpos{k};
                rawtracking(k).data(:,end+1)=k; % tack on the epoch to the tracking data
            end
            thisDay.coords=cell2mat({rawtracking.data}');
            % now we can register the tracking with thisday.epochs to see
            % when it is
            fprintf('extracted tracking from %s, %d runs, duration: %.2f min \n',...
                mousefiles(msindex).tracking{posfile},length(unique(thisDay.coords(:,end))),thisDay.coords(end,1)/60);
            
        end
        
        posfile=find(contains(mousefiles(msindex).tracking,['pos' daynum]) & ...
            contains(mousefiles(msindex).tracking,'lin'));
        if ~isempty(posfile)
            rawtracking=load(fullfile(mousefiles(msindex).directory,mousefiles(msindex).tracking{posfile}));
            thisDay.rawlintracking=rawtracking.linpos(~isempty(rawtracking.linpos));
            thisDay.sourcefiles=[thisDay.sourcefiles; {fullfile(mousefiles(msindex).directory,mousefiles(msindex).tracking{posfile})}];
            epochpos=rawtracking.linpos{cellfun(@(a) ~isempty(a),rawtracking.linpos)};
            
            clear rawtracking;
            for k=1:length(epochpos)
                rawtracking(k)=epochpos{k};
                rawtracking(k).data(:,end+1)=k; % tack on the epoch to the tracking data
            end
            thisDay.lincoords=cell2mat({rawtracking.data}');
            % now we can register the tracking with thisday.epochs to see
            % when it is
            fprintf('extracted tracking from %s, %d runs, duration: %.2f min \n',...
                mousefiles(msindex).tracking{posfile},length(unique(thisDay.coords(:,end))),thisDay.coords(end,1)/60);
        end
        
        
        
        
        % linearized position looks like it hasnt been caluclated for very
        % many sessions, so im just going to ignore it for now. I'll likely
        % go through each session and given its type and track, draw the
        % linear tracks myself, however i'll probably have to create a new
        % way to track for the alternation maze
        %{
         linposfile=find(contains(mousefiles(i).tracking,['pos' daynum]) & ...
            contains(mousefiles(i).tracking,'lin'));
        if ~isempty(linposfile)
            rawlinear=load(fullfile(mousefiles(i).directory,mousefiles(i).tracking{linposfile}));
            rawlin=rawlinear.linpos(~isempty(rawlinear.linpos));
            epochpos=rawlin{1};
            clear rawlinear;
            for k=1:length(epochpos)
                rawlinear(k)=epochpos{k};
                rawlinear(k).data(:,end+1)=k; % tack on the epoch to the tracking data
            end
            thisDay.rawtracking.linear=rawlinear;
        end
        %}
        % find this days spikefile
        spikefile=find(contains(mousefiles(msindex).unitdata,['spikes' daynum]) );
        if ~isempty(spikefile)
            try
                rawunits=load(fullfile(mousefiles(msindex).directory,mousefiles(msindex).unitdata{spikefile}));
                thisDay.sourcefiles=[thisDay.sourcefiles; {fullfile(mousefiles(msindex).directory,mousefiles(msindex).unitdata{spikefile})}];
                % get in to the epoch cells
                rawunits.spikes=rawunits.spikes{cellfun(@(a) ~isempty(a), rawunits.spikes)};
                
                thisDay.oldspikes=rawunits;
                
                
                % now go in and grab the unit data (across small epochs)...
                cumct=1; clear unitdata;
                tempstruct=struct('tet',i,'unitnum',j,'ts',[],'meanrate',[],'tag',[],'fields',[],'descript',[],'epoch',[]);
                
                % this just concatenates all the cells so we can mush them
                % together later
                for i=1:length(rawunits.spikes) % always full
                    % for each tetrode
                    for j=1:length(rawunits.spikes{i}) % sometimes empty
                        if ~isempty(rawunits.spikes{i}{j})
                            % for unit on that tetrode
                            for k=1:length(rawunits.spikes{i}{j}) % often empty
                                if ~isempty(rawunits.spikes{i}{j}{k})
                                    %fprintf('index: epoch %d, tetrode %d, unit %d  had mean rate %.2f \n',i, j, k, rawunits.spikes{i}{j}{k}.meanrate);
                                    unitdata(cumct)=tempstruct; % load the blank struct
                                    unitdata(cumct).ts=rawunits.spikes{i}{j}{k}.data; % add data
                                    unitdata(cumct).ts(:,end+1)=i; % add which epoch that spike is from
                                    unitdata(cumct).tet=j; unitdata(cumct).unitnum=k; % tet and unit
                                    unitdata(cumct).tag= rawunits.spikes{i}{j}{k}.tag; % the tag
                                    unitdata(cumct).meanrate=rawunits.spikes{i}{j}{k}.meanrate; % mean rate
                                    unitdata(cumct).fields=rawunits.spikes{i}{j}{k}.fields; % fields
                                    unitdata(cumct).descript=rawunits.spikes{i}{j}{k}.descript; % and descript
                                    unitdata(cumct).epoch=i;
                                    cumct=cumct+1;
                                end
                            end
                        else
                            %fprintf('sess %d, tetrode %d was empty \n',i,j)
                        end
                    end
                end
                fprintf('extracted Units from %s \n',mousefiles(msindex).unitdata{spikefile});
                
                % now we collapse them into single units
                if exist('unitdata','var')
                    activetets=unique([unitdata.tet]);
                    clear newunitdata; cumct=1;
                    for i=1:length(activetets)
                        cellcts=unique([unitdata([unitdata.tet]==activetets(i)).unitnum]);
                        for j=1:length(cellcts)
                            newunitdata(cumct)=unitdata(find([unitdata.tet]==activetets(i) & [unitdata.unitnum]==cellcts(j),1,'first'));
                            newunitdata(cumct).ts=cell2mat({unitdata([unitdata.tet]==activetets(i) & [unitdata.unitnum]==cellcts(j)).ts}');
                            cumct=cumct+1;
                        end
                    end
                end
                % animdb recording data.tetinfo to get the tetrode
                % locations
                
                
                % if things didnt work...
                if exist('newunitdata','var')
                    thisDay.unitdata=newunitdata;
                    fprintf('Organized Units from %s, there are %d \n',mousefiles(msindex).unitdata{spikefile}, length(newunitdata));
                else
                    keyboard
                    thisDay.unitdata=[];
                    fprintf('There were no units from %s \n',mousefiles(msindex).unitdata{spikefile});
                end
                % if things really didnt work...
            catch
                fprintf('couldnt extract Units from %s \n',mousefiles(msindex).unitdata{spikefile});
                keyboard
            end
        else
            keyboard
        end
        
        %%%%%%%%%%% now at least find the LFP files so we can grab them
        %%%%%%%%%%% when we want.. we'll prob grab the tet with the most
        %%%%%%%%%%% units on it
        
        LFPfiles=mousefiles(msindex).eeg;
        sessLFP=LFPfiles(contains(LFPfiles,['eeg' daynum]));
        thisDay.LFPfiles=sessLFP;
        % I'll analyze my own LFP data...
        
        
        
        
        SuperMouse(sesstally)=thisDay; sesstally=sesstally+1;
        fprintf('added %s %s to supermouse \n', thisDay.name, thisDay.day);
    end
end


%%
% okay lets see if we can at least peek at the tracking data
for i=1:10
    if ~isempty(SuperMouse(i).coords)
        % parse the epochs
        figure;
        [epochs,~,epinds]=unique({SuperMouse(i).epochs.epoch_type});
        for k=1:length(epochs)
            % grab epoch indices that arent this epoch
            nanepochs=[SuperMouse(i).epochs(epinds~=k).epoch];
            
            subplot(1,length(epochs),k);
            tempcoords=SuperMouse(i).coords;
            for no=1:length(nanepochs), tempcoords(tempcoords(:,7)==nanepochs(no),:)=nan; end
            plot(tempcoords(:,2),tempcoords(:,3));
            title(epochs{k});
        end
        sgtitle(sprintf('Mouse %s day %s',SuperMouse(i).name, SuperMouse(i).daynum));
    end
    
end
%%
for i=10
    if ~isempty(SuperMouse(i).coords)
        % parse the epochs
        %figure;
        [epochs,~,epinds]=unique({SuperMouse(i).epochs.epoch_type});
        for k=1:length(epochs)
            % grab epoch indices that arent this epoch
            nanepochs=[SuperMouse(i).epochs(epinds~=k).epoch];
            tempcoords=SuperMouse(i).coords;
            for no=1:length(nanepochs), tempcoords(tempcoords(:,7)==nanepochs(no),:)=nan; end
            % now use that coordinate data to plot a cells firing fields
            tempcoords=tempcoords(~isnan(tempcoords(:,1)),1:3);
            
            for un=1:length(SuperMouse(i).unitdata)
                Session=struct('edit_coords',tempcoords);
                unit=struct('ts',SuperMouse(i).unitdata(un).ts(:,1));
                [a,b,c]=cell_SmoothPlacePlot(Session,unit,'Factor',2);
            end
            sgtitle(sprintf('Mouse %s day %s',SuperMouse(i).name, SuperMouse(i).daynum));
            
        end
    end
    
end

%% this is because I want velocity data and its a bitch to get
velsmooth=.2; PixPerCM=6; timejumpdelay=1;
for i=62:length(SuperMouse)
    tic
    if ~isempty(SuperMouse(i).coords)
        % parse the epochs
        %figure;
        [epochs,~,epinds]=unique({SuperMouse(i).epochs.epoch_type});
        for k=1:length(epochs)
            if contains(epochs(k),'Track')
                % grab epoch indices that arent this epoch
                nanepochs=[SuperMouse(i).epochs(epinds~=k).epoch];
                tempcoords=SuperMouse(i).coords;
                for no=1:length(nanepochs), tempcoords(tempcoords(:,7)==nanepochs(no),:)=nan; end
                
                % now get our velocity data
                tempcoords=tempcoords(~isnan(tempcoords(:,1)),1:3);
                [~,uniquets]=unique(tempcoords(:,1));
                tempcoords=tempcoords(uniquets,:);
                
                dt = diff(tempcoords(:,1)); % difference between consecutive coords
                % chop off the first timestamp because we dont know its elapsed time
                timestamps = tempcoords((2:end),1); % remove coords so dt1 is for before tcoord1
                xycoords = tempcoords((2:end),2:3); % save raw xy data (because we smooth for velocity filtering)
                
                
                
                % this is the smoothing algorithm to get the velocity data
                windowspan=round(velsmooth/median(dt))*5; % how many seconds to smooth the data across?
                coorddata2=[SmoothMat2(tempcoords(:,2),[1 windowspan*4], windowspan),...
                    SmoothMat2(tempcoords(:,3),[1 windowspan*4], windowspan)];
                
                
                % this is in cm
                displacement=sqrt(diff(coorddata2(:,1)).^2+diff(coorddata2(:,2)).^2)/PixPerCM;
                %this is a non smoothed displacement, we may want to bin this to remove the
                %small jitters that occur
                % any displacement where the time is too far, we have to nan out
                elapsed=dt;
                % nan out long lags so they dont contribute to the smoothed average
                displacement(elapsed>timejumpdelay)=nan;
                elapsed(elapsed>timejumpdelay)=nan;
                % for smoothing the velocity data, we grab all the coords that are within
                % our timewindow average them and thats our smoothed velocity.
                velocity=nan(length(displacement),1);
                
                % this takes fucking forever, so lets just do a 5 bin arma
                parfor q=1:length(velocity)
                    velocity(q)=nansum(displacement(abs(timestamps-timestamps(q))<=velsmooth))/...
                        nansum(elapsed(abs(timestamps-timestamps(q))<=velsmooth));
                    
                end
                SuperMouse(i).TrackVelData=velocity;
                fprintf('Session %d %s %s done in %.2f mins \n', i,...
                    SuperMouse(i).name, SuperMouse(i).daynum, toc/60);
                
                
            end
        end
    end
end









%% Clean up workspace...
clearvars -except SuperMouse mousefiles