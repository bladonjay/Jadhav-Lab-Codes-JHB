
clear

topDir = cs_setPaths();

%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animals = {'CS33'};
minstd = 3;

%change minnum depending on whether ripples are found using tetrodes that
%have cells, or single tetrode with the most cells
minnum = 1; 
velfilter = 4; %velthresh
%%
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt',filesep,animal,'_direct',filesep];
    
    tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
    runepochs = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(runepochs(:,1));
    
    for d = 1:length(days)
        day = days(d);
        dayepochs = runepochs(runepochs(:,1) == day,:);
        
        %get ripple times
        %clear rippletime
        rippletime = DFTFsj_getriptimes(animDir,animal, dayepochs,'tetfilter', '(isequal($descrip, ''riptet''))','minthresh',minstd);
        
        %% find periods of immobility
        if ~isempty(velfilter)
            pos = loaddatastruct(animDir, animal, 'pos', day); % get animal [time, positions, velocity]
        end
        epochs = dayepochs(:,2)';
        for ep = epochs
            if isempty(rippletime)
                rip.timevec = [];
                rip.starttime =[];
                rip.endtime = [];
                rip.ripvec = [];
                rip.total_duration = 0 ;
                
                rip.minthresh = minstd;
                rip.minnumber = minnum;
                rip.velthresh = velfilter;
                ripple{day}{ep} = rip;
                continue
            end
            %ripple periods
            riplist = vec2list(rippletime{day}{ep}.nripples >= minnum,rippletime{day}{ep}.time);
            if ~isempty(velfilter)
                if size(pos{day}{ep}.data,2) > 5  % get velocity
                    velocity = pos{day}{ep}.data(:,9);% smoothed velocity
                else
                    velocity = pos{day}{ep}.data(:,5);
                end
                postimes = pos{day}{ep}.data(:,1);  % get time stamps
                immobile = vec2list(velocity < velfilter,postimes); % generate [start end] list of immobile epochs
            end
            
            if ~isempty(riplist)
                times_filteeg = rippletime{day}{ep}.time;
                duration = (times_filteeg(end) - times_filteeg(1)).*1000; % ms
                
                [imtime,im_vec] = wb_list2vec(immobile,times_filteeg);
                [riptime,ripvec] = wb_list2vec(riplist,times_filteeg);
                nonripvec = ~ripvec;
                
                
                if ~isempty(velfilter)
                    ripvec = ripvec & im_vec;
                    riplist = vec2list(ripvec,times_filteeg);
                    nonripvec = ~ripvec & ~im_vec;
                end
                nriplist = vec2list(nonripvec,times_filteeg);
                rip.timevec = times_filteeg;
                rip.starttime = riplist(:,1);
                rip.endtime = riplist(:,2);
                rip.ripvec = ripvec;
                rip.nonripvec = nonripvec;
                rip.nstarttime =  nriplist(:,1);
                rip.nendtime =  nriplist(:,2);
                ripduration = round(sum(riplist(:,2)-riplist(:,1)));
                rip.total_duration = ripduration;
                
                %                 im.timevec = times_filteeg;
                %                 im.starttime = immobile(:,1);
                %                 im.endtime = immobile(:,2);
                %                 im.total_duration =  round(sum(immobile(:,2)-immobile(:,1)));
            else
                rip.timevec = [];
                rip.starttime =[];
                rip.endtime = [];
                rip.ripvec = [];
                rip.total_duration = 0 ;
            end
            rip.minthresh = minstd;
            rip.minnumber = minnum;
            rip.velthresh = velfilter;
            ripple{day}{ep} = rip;
            
            %             im.velthresh = velfilter;
            %             immobility{day}{ep} = im;
            sprintf('a %s d %d e %d: %d s of ripples',animal,day,ep,ripduration)
            %             sprintf('d %d e %d: %d s of immobility',day,ep,round(sum(immobile(:,2)-immobile(:,1))))
            clear rip;clear im;
        end
    end
    save(sprintf('%s/%srippletimes.mat', animDir, animal), 'ripple');
    %save(sprintf('%s/%simmobiletime%02d.mat', animaldir, animalprefix, day), 'immobility');
    clear ripple;clear immobility
    
end
