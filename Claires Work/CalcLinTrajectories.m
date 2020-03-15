%% CalcLinTrajectories

% this gathers the linear place fields for each cell, and calcs their
% information score and field parameters (field size, peak rate, and
% sparsity



% for it to be a place field we'll use standard thresholds:
% 2. the peak must be above 2 hz
% 3. the field (top 75% of field using flood fill) must be less than 50% of
%    the maze  its on line 279


%%  first we gather the linearized trajectories for each unit
% we want to do this just for the outbound trajectories, and it can be a
% place field for either (all the cells will have a left and a right out
% field

savedir=uigetdir;
runboots=250;
verbose=1; % start with true, move to false
runblocks=0;

velthreshold=3;
velsmooth=8;
% gonna add a few fields to the units struct
FieldPropsNan=struct('PFmax',nan(1,4),'PFmaxpos',nan(1,4),...
    'info',nan(1,4),'sparsity',nan(1,4),...
    'Zpfmax',nan(1,4),'PFsize',nan(1,4),...
    'PFmaxP',nan(1,4),'infoP',nan(1,4),...
    'sparsityP',nan(1,4));

% and also


for ses=1:length(SuperRat)
    tic
    if SuperRat(ses).longTrack
        % 1 is bottom left, 2 is bottom right, 3 is home
        trajinds=[3 1; 3 2; 1 3; 2 3];
        % 45 are the left linearized positions
        % 67 are the right linearized positions
        keepinds=[];
        allposplot=[];
        cellnotes={SuperRat(ses).units.tag};
        for cn=1:length(cellnotes)
            if isempty(cellnotes{cn})
                cellnotes(cn)={'mua'};
            end
        end
        % grab the run blocks
        blocks=unique(SuperRat(ses).LinCoords(:,6));
        fulltraj=SuperRat(ses).LinCoords;
        if isfield(SuperRat(ses).units(1),'LinSpikeCts')
            SuperRat(ses).units=rmfield(SuperRat(ses).units,{'LinSpikeCts',...
                'LinPlaceFields','FiresDuringRun','PFexist','RunRates'});
        end
        % run vs sleep blocks
        for j=1:length(SuperRat(ses).units)
            
            
            spikes=SuperRat(ses).units(j).ts;
            try
                SuperRat(ses).units(j)=rmfields(SuperRat(ses).units(j),...
                    {'LinPlaceFields','LinSpikeCts','FiresDuringRun','PFexist','RunRates'});
            end
            
            SuperRat(ses).units(j).FieldProps=FieldPropsNan; % initialize the field props struct
            if runblocks
                for bl=1:length(blocks)
                    % for each trajectory
                    for tr=1:4
                        % this block and this trajectory
                        keepinds=fulltraj(:,8)==trajinds(tr,1) & fulltraj(:,9)==trajinds(tr,2) & fulltraj(:,10)==blocks(bl);
                        % are there any indices to match? if not nan everything out
                        if sum(keepinds)>30*15 % 15 seconds of data
                            temptraj=fulltraj(keepinds,[1:3 linpull(tr,:)]); % pull fav line, and rows
                            
                            % and nan out the slow times (dx for linear pos)
                            temptraj=sortrows(temptraj,1); % order chronologically
                            tooslow=SmoothMat2(abs(diff(temptraj(:,5))),[0 50],10)<1;
                            temptraj(tooslow,:)=[];
                            
                            % find the epochs so you can kill bad spikes
                            breaks=find(diff(temptraj(:,1))>1);
                            epochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1); temptraj(end,1)]];
                            % only use epochs that last more than a seconds
                            epochs(diff(epochs,1,2)<=1,:)=[];
                            % get occupancy
                            [bincts,bins]=histcounts(temptraj(:,5),100); % in bins
                            occupancy=bincts/30;
                            smoothoccup=SmoothMat2(occupancy,[5 0],2);
                            SuperRat(ses).LinOccup(tr,:)=smoothoccup;
                            fullcurves={}; cellsort=[];
                            % now for each unit capture a mean tuning curve for each run
                            try
                                Espikes=EpochCoords(spikes(:,1),epochs); % only pull spike ts in the run epochs
                                [~,Erates]=event_spikes(spikes,epochs(:,1),0,epochs(:,2));
                                spikepos=interp1(temptraj(:,1),temptraj(:,5),Espikes,'nearest'); % now interp to position
                                %if contains(cellnotes{j},'accepted')
                                % need to bootstrap here so i can get a p value on peak
                                % width and information score (also could do positional
                                % info for each pass...
                                if verbose, figure; end
                                
                                if runboots>0
                                    ha=waitbar(0,'Starting');
                                    if~isempty(Espikes)
                                        Bmax=nan(1,runboots); Binfo=nan(1,runboots); Bsparsity=nan(1,runboots);
                                        for bt=1:runboots
                                            EBspikes=[];
                                            for i=1:length(epochs)
                                                thesespikes=EpochCoords(spikes(:,1),epochs(i,:));
                                                % circ shift the spikes in this epoch by a min of 1/2 second
                                                EBspikes=[EBspikes; PermuteSpikeTimes(thesespikes,[epochs(i,1); epochs(i,2)],0.5)];
                                            end
                                            Bspikepos=interp1(temptraj(:,1),temptraj(:,5),EBspikes,'nearest');
                                            bnspikes=histcounts(Bspikepos,bins); % get number of spikes per position
                                            smoothspikes=SmoothMat2(bnspikes,[5 0],2); % two pixel kernel
                                            % calculate null info
                                            [Binfo(bt),Bsparsity(bt)]=Skaggs_basic(smoothoccup./max(smoothoccup),...
                                                smoothspikes,nanmean(smoothspikes));
                                            % calculate null peak rate
                                            bLinPlaceField=SmoothMat2(bnspikes./occupancy,[5 0],2);
                                            if verbose, plot(bLinPlaceField); hold on; end
                                            Bmax(bt)=max(bLinPlaceField);
                                            waitbar(bt/runboots,ha,'working on it ...');
                                        end
                                    elseif isempty(Espikes)
                                        Binfo=nan; Bsparsity=nan; Bmax=nan;
                                    end
                                    close(ha);
                                end
                                
                                [nspikes,bins]=histcounts(spikepos,bins); % get number of spikes per position
                                smoothspikes=SmoothMat2(nspikes,[5 0],2);
                                % now save this tuning curve
                                LinPlaceField=SmoothMat2(nspikes./occupancy,[5 0],2);
                                if verbose, plot(LinPlaceField,'LineWidth',3); end
                                
                                % now qualify this as a pf or not
                                [pfmax,pfmaxpos]=max(LinPlaceField); % has to be above 5
                                Zpfmax=max(zscore(LinPlaceField)); % has to be above 3
                                
                                % run skaggs on individually smoothed occupancy and
                                % spikes
                                [info,sparsity]=Skaggs_basic(smoothoccup./max(smoothoccup),...
                                    smoothspikes,nanmean(smoothspikes));
                                
                                % and now pf size e.g. size of top 75%% of max rate
                                pfstart=find(LinPlaceField(1:pfmaxpos) < pfmax*.25,1,'last');
                                if isempty(pfstart), pfstart=1; end
                                pfend=pfmaxpos+find(LinPlaceField(pfmaxpos:end) < pfmax*.25,1,'first');
                                if isempty(pfend), pfend=100; end
                                
                                % try to calculate pf size, otherwise it doesnt exist
                                pfsize=pfend-pfstart;
                                
                                
                                % now add to the session struct
                                
                                
                                %figure; plot(LinPlaceField); hold on;
                                %plot([pfstart pfstart],[0 LinPlaceField(pfstart)],'r');
                                %plot([pfend pfend],[0 LinPlaceField(pfend)],'r');
                                SuperRat(ses).units(j).LinPlaceFields{bl}(tr,:)=LinPlaceField;
                                SuperRat(ses).units(j).LinSpikeCts{bl}(tr,:)=smoothspikes;
                                SuperRat(ses).units(j).FiresDuringRun(bl,tr)=length(Espikes)>size(epochs,1); % is this cell active?
                                SuperRat(ses).units(j).PFexist(bl,tr)=pfmax>3 && pfsize<50; % is it a place cell
                                SuperRat(ses).units(j).RunRates{bl,tr}=Erates; % run the splitter score later
                                
                                SuperRat(ses).units(j).FieldProps.PFmax(bl,tr)=pfmax; % peak rate
                                SuperRat(ses).units(j).FieldProps.PFmaxpos(bl,tr)=pfmaxpos; % position on maze
                                SuperRat(ses).units(j).FieldProps.info(bl,tr)=info; % info
                                SuperRat(ses).units(j).FieldProps.sparsity(bl,tr)=sparsity;
                                SuperRat(ses).units(j).FieldProps.Zpfmax(bl,tr)=Zpfmax; % zscored max rate
                                SuperRat(ses).units(j).FieldProps.PFsize(bl,tr)=pfsize; % size in % of trajectory
                                if runboots>0
                                    SuperRat(ses).units(j).FieldProps.PFmaxP(bl,tr)=1-normcdf(pfmax,nanmean(Bmax),nanstd(Bmax)); % size in % of trajectory
                                    SuperRat(ses).units(j).FieldProps.infoP(bl,tr)=1-normcdf(info,nanmean(Binfo),nanstd(Binfo)); % size in % of trajectory
                                    SuperRat(ses).units(j).FieldProps.sparsityP(bl,tr)=1-normcdf(sparsity,nanmean(Bsparsity),nanstd(Bsparsity)); % size in % of trajectory
                                end
                                if verbose, title(sprintf('Info: %.2f P=%.5f',info,normpdf(info,nanmean(Binfo),nanstd(Binfo)))); end
                                
                                
                                % end
                            catch
                                fprintf('Sess %d block %d Trajectory %d cell %d couldnt be analyzed \n',ses, bl, tr, j);
                            end % end of try
                        end % if we have trajectories
                    end % of trajectories
                end % blocks
                %  now one for all the blocks
            else % if run whole day as unitary session
                for tr=1:4
                    % this this trajectory (e.g. start and finish are
                    % correct
                    keepinds=fulltraj(:,4)==trajinds(tr,1) & fulltraj(:,5)==trajinds(tr,2);
                    % did the animal spend enough time in this trajectory?
                    if sum(keepinds)>30*15 % 15 seconds of data
                        temptraj=fulltraj(keepinds,:); % pull fav line, and rows
                        
                        % and nan out the slow times (dx for linear pos)
                        % thjis is the velocity threshold
                        temptraj=sortrows(temptraj,1); % order chronologically
                        % 50 bin span, but a std of only 1/4 second
                        % or like 3/4 of a second
                        tooslow=SmoothMat2(temptraj(:,7),[0 50],velsmooth)<=velthreshold;
                        temptraj(tooslow,:)=[];
                        
                        % find the epochs so you can kill bad spikes
                        breaks=find(diff(temptraj(:,1))>1);
                        epochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1); temptraj(end,1)]];
                        % only use epochs that last more than a seconds
                        epochs(diff(epochs,1,2)<=1,:)=[];
                        % get occupancy
                        [bincts,bins]=histcounts(temptraj(:,8),1:100); % in bins
                        occupancy=bincts/30;
                        smoothoccup=SmoothMat2(occupancy,[5 0],2);
                        SuperRat(ses).LinOccup(tr,:)=smoothoccup;
                        fullcurves={}; cellsort=[];
                        % now for each unit capture a mean tuning curve for each run
                        try
                            Espikes=EpochCoords(spikes(:,1),epochs); % only pull spike ts in the run epochs
                            [~,Erates]=event_spikes(spikes,epochs(:,1),0,epochs(:,2));
                            spikepos=interp1(temptraj(:,1),temptraj(:,8),Espikes,'nearest'); % now interp to position
                            %if contains(cellnotes{j},'accepted')
                            % need to bootstrap here so i can get a p value on peak
                            % width and information score (also could do positional
                            % info for each pass...
                            if verbose, figure; end
                            
                            if runboots>0
                                %ha=waitbar(0,'Starting');
                                if~isempty(Espikes)
                                    Bmax=nan(1,runboots); Binfo=nan(1,runboots); Bsparsity=nan(1,runboots);
                                    for bt=1:runboots
                                        EBspikes=[];
                                        for i=1:length(epochs)
                                            % circ shift the spikes in this epoch by a min of 1/2 second
                                            EBspikes=[EBspikes; PermuteSpikeTimes(EpochCoords(spikes(:,1),epochs(i,:)),[epochs(i,1); epochs(i,2)],0.25)];
                                        end
                                        Bspikepos=interp1(temptraj(:,1),temptraj(:,8),EBspikes,'nearest');
                                        bnspikes=histcounts(Bspikepos,bins); % get number of spikes per position
                                        smoothspikes=SmoothMat2(bnspikes,[5 0],2); % two pixel kernel
                                        % calculate null info
                                        [Binfo(bt),Bsparsity(bt)]=Skaggs_basic(smoothoccup./max(smoothoccup),...
                                            smoothspikes,nanmean(smoothspikes));
                                        % calculate null peak rate
                                        bLinPlaceField=SmoothMat2(bnspikes./occupancy,[5 0],2);
                                        if verbose, plot(bLinPlaceField); hold on; end
                                        Bmax(bt)=max(bLinPlaceField);
                                        %waitbar(bt/runboots,ha,'working on it ...');
                                    end
                                elseif isempty(Espikes)
                                    Binfo=nan; Bsparsity=nan; Bmax=nan;
                                end
                                %close(ha);
                            end
                            
                            [nspikes,bins]=histcounts(spikepos,bins); % get number of spikes per position
                            smoothspikes=SmoothMat2(nspikes,[5 0],2);
                            % now save this tuning curve
                            % smooth the rate map AFTER dividing by raw
                            % occupancy
                            LinPlaceField=SmoothMat2(nspikes./occupancy,[5 0],2);
                            if verbose
                                plot(LinPlaceField,'LineWidth',3);
                            end
                            
                            % now qualify this as a pf or not
                            [pfmax,pfmaxpos]=max(LinPlaceField); % has to be above ??
                            Zpfmax=max(nanzscore(LinPlaceField)); % has to be above ??
                            
                            % run skaggs on individually smoothed occupancy and
                            % spikes
                            [info,sparsity]=Skaggs_basic(smoothoccup./max(smoothoccup),...
                                smoothspikes,nanmean(smoothspikes));
                            
                            % and now pf size e.g. size of top 75%% of max rate
                            % this will underestimate, because it finds the
                            % first that drops below the high rate 'flood
                            % fill' method
                            pfstart=find(LinPlaceField(1:pfmaxpos) < pfmax*.25,1,'last');
                            if isempty(pfstart), pfstart=1; end
                            pfend=pfmaxpos+find(LinPlaceField(pfmaxpos:end) < pfmax*.25,1,'first');
                            if isempty(pfend), pfend=100; end
                            
                            % try to calculate pf size, otherwise it doesnt exist
                            pfsize=pfend-pfstart;
                            
                            % now add to the session struct
                            SuperRat(ses).units(j).LinPlaceFields{1}(tr,:)=LinPlaceField;
                            SuperRat(ses).units(j).LinSpikeCts{1}(tr,:)=smoothspikes;
                            SuperRat(ses).units(j).FiresDuringRun(1,tr)=pfmax>1; % is this cell active? e.g. does it pass 1 hz at anywhere on track
                            SuperRat(ses).units(j).RunRates{1,tr}=Erates; % run the splitter score later
                            
                            SuperRat(ses).units(j).FieldProps.PFmax(1,tr)=pfmax; % peak rate
                            SuperRat(ses).units(j).FieldProps.PFmaxpos(1,tr)=pfmaxpos; % position on maze
                            SuperRat(ses).units(j).FieldProps.info(1,tr)=info; % info
                            SuperRat(ses).units(j).FieldProps.sparsity(1,tr)=sparsity;
                            SuperRat(ses).units(j).FieldProps.Zpfmax(1,tr)=Zpfmax; % zscored max rate
                            SuperRat(ses).units(j).FieldProps.PFsize(1,tr)=pfsize; % size in % of trajectory
                            if runboots>0
                                SuperRat(ses).units(j).FieldProps.PFmaxP(1,tr)=1-normcdf(pfmax,nanmean(Bmax),nanstd(Bmax)); % peak statistically high?
                                SuperRat(ses).units(j).FieldProps.infoP(1,tr)=1-normcdf(info,nanmean(Binfo),nanstd(Binfo)); % information statistically high?
                                SuperRat(ses).units(j).FieldProps.sparsityP(1,tr)=1-normcdf(sparsity,nanmean(Bsparsity),nanstd(Bsparsity)); % Sparsity?
                                SuperRat(ses).units(j).PFexist(1,tr)=1-normcdf(pfmax,nanmean(Bmax),nanstd(Bmax))<.05 && Zpfmax>2 && pfsize<75; % is it a place cell
                            else
                                SuperRat(ses).units(j).PFexist(1,tr)=PFmax>2 && Zpfmax>2 && pfsize<75; % is it a place cell
                            end

                            if verbose
                                title(sprintf('Info: %.2f P=%.5f Peak: %.2f, p=%.5f',info,1-normcdf(info,nanmean(Binfo),nanstd(Binfo)),...
                                    pfmax,1-normcdf(pfmax,nanmean(Bmax),nanstd(Bmax)))); 
                                answer = questdlg('Stop verbose?', 'Verbose Menu', 'Yes, quiet','no, more figs','Yes, quiet');
                                verbose=contains(answer,'no');
                            end
                        catch
                            fprintf('Whole Sess %d Trajectory %d cell %d couldnt be analyzed \n',ses, tr, j);
                        end % end of try
                    end % if we have trajectories
                end % of trajectories
            end
            fprintf('%d ',j);
        end % units
        fprintf('\n Session %d done in %.2f mins \n',ses,toc/60);
    end
end
if any(savedir~=0)
    save(fullfile(savedir,'SuperRatAllPFStats'),'SuperRat');
    fprintf('Finished Boot and Saved out \n');
else
    fprintf('finished boot, didnt save \n');
end
%% this runs the comparison of mean rates for each trajectory
% Basically its a low-res splitter score, if it doesnt work we'll have to
% do a spatial cross correlation between the trajectories


boot=500;
for i=1:length(SuperRat)
    tic
    ha=waitbar(0,'Starting Session');
    try SuperRat(i).units=rmfield(SuperRat(i).units,'TrajScores'); end
    for j=1:length(SuperRat(i).units)
        
        for bl=1:size(SuperRat(i).units(j).RunRates,1)
            RunCell=SuperRat(i).units(j).RunRates(bl,:);
            RunRates=[RunCell{1} RunCell{2} RunCell{3} RunCell{4};...
                ones(1,length(RunCell{1})) ones(1,length(RunCell{2}))+1 ...
                ones(1,length(RunCell{3}))+2 ones(1,length(RunCell{4}))+3]';
            % out left vs out right (1, 2), in left vs in right
            % positive is left preferring, negative is right preferring
            % for now these scores are total spikes, maybe well need max
            % which would take a different math here for the P value
            TrajScore=(nanmean(RunRates(RunRates(:,2)==1,1))-nanmean(RunRates(RunRates(:,2)==2,1)))/...
                (nanmean(RunRates(RunRates(:,2)==1,1))+nanmean(RunRates(RunRates(:,2)==1,1)));
            TrajScore(3)=(nanmean(RunRates(RunRates(:,2)==3,1))-nanmean(RunRates(RunRates(:,2)==4,1)))/...
                (nanmean(RunRates(RunRates(:,2)==3,1))+nanmean(RunRates(RunRates(:,2)==4,1)));
            
            if boot>0
                % get the pull for
                outTraj=RunRates(RunRates(:,2)==1 | RunRates(:,2)==2,1);
                inTraj=RunRates(RunRates(:,2)==3 | RunRates(:,2)==4,1);
                TrajBoot=[nan nan];
                for bt=1:boot
                    % randomize the order so you can even odd them
                    outOrder=randperm(length(outTraj));
                    inOrder=randperm(length(inTraj));
                    TrajBoot(bt,1)=abs(nanmean(outTraj(mod(outOrder,2)==1))-...
                        nanmean(outTraj(mod(outOrder,2)==0)))/...
                        (nanmean(outTraj(mod(outOrder,2)==1))+...
                        nanmean(outTraj(mod(outOrder,2)==0)));
                    TrajBoot(bt,2)=abs(nanmean(inTraj(mod(inOrder,2)==1))-...
                        nanmean(inTraj(mod(inOrder,2)==0)))/...
                        (nanmean(inTraj(mod(inOrder,2)==1))+...
                        nanmean(inTraj(mod(inOrder,2)==0)));
                end
                TrajScore(2)=1-normcdf(abs(TrajScore(1)), nanmean(TrajBoot(:,1)),nanstd(TrajBoot(:,1)));
                TrajScore(4)=1-normcdf(abs(TrajScore(3)), nanmean(TrajBoot(:,2)),nanstd(TrajBoot(:,2)));
            else
                TrajScore(2)=nan; TrajScore(4)=nan;
            end
            SuperRat(i).units(j).TrajScores(bl,:)=TrajScore;
        end
        waitbar(j/length(SuperRat(i).units),ha,sprintf('running unit %d',j));
    end
    close(ha);
    fprintf('Session %d took %.2f minutes \n',i,toc/60);
    
end


%% just look at some cells to see if they're selective and that
% the raw rates look reasonable


% now lets see if there are valid priors here
for i=1:length(SuperRat(ses).units)
    if any(SuperRat(ses).units(i).PFexist)
        figure;
        for j=1:length(SuperRat(ses).RunEpochs)
            subplot(1,length(SuperRat(ses).RunEpochs),j)
            plot(SuperRat(ses).units(i).LinPlaceFields{j}(1,:));
            hold on;
            plot(-SuperRat(ses).units(i).LinPlaceFields{j}(2,:));
            plot(SuperRat(ses).units(i).LinPlaceFields{j}(3,:));
            plot(-SuperRat(ses).units(i).LinPlaceFields{j}(4,:));
            
            title([num2str(i) ' ' SuperRat(ses).units(i).area SuperRat(ses).units(i).type]);
        end
    end
end
% looks really good actually!

%%


