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
runboots=500;
verbose=0; % start with true, move to false
runblocks=0;
calctheta=false; % do you want to calculate theta pref and precession slope?

trajinds=[3 1; 3 2; 1 3; 2 3];

speedthreshold=3;
veltimesmooth=8; % bins
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
        % 4,5 are the left linearized positions
        % 6,7 are the right linearized positions
        keepinds=[]; allposplot=[];
        % grab the run blocks
        blocks=SuperRat(ses).RunEpochs;
        fulltraj=SuperRat(ses).LinCoords;
        fulltraj(~ismember(fulltraj(:,6),blocks),:)=[]; % use only run sesses
        epochbreaks=fulltraj(diff(fulltraj(:,6))~=0,1);
        if any(epochbreaks)
            epochtimes=[fulltraj(1) epochbreaks'+1; epochbreaks' fulltraj(length(fulltraj),1)]';
        else
            epochtimes=[fulltraj(1) fulltraj(end,1)];
        end
        
        if calctheta
            thetadata=load(fullfile(SuperRat(ses).LFP.filedir,SuperRat(ses).LFP.filename));
            thetadata=thetadata.eegstruct; eegmat=[];
            for k=1:length(thetadata)
                eegmat=[eegmat; double(thetadata(k).starttime:1/thetadata(i).filtersamprate:thetadata(k).endtime)'...
                    double(thetadata(k).data)];
            end
            eegmat(:,3)=rescale(eegmat(:,3),-pi,pi); clear thetadata;
        end

        % run blocks
        for j=1:length(SuperRat(ses).units)
            FieldProps=FieldPropsNan; % initialize the field props struct
            spikes=SuperRat(ses).units(j).ts;
            [allspikes,epochspk]=event_spikes(spikes,epochtimes(:,1),0,diff(epochtimes,1,2));           

            % remove all runs in epochs with no spikes
            keepblocks=blocks(epochspk==0);
            celltraj=fulltraj(ismember(fulltraj(:,6),keepblocks),:);
            clear thetastats;
            
            % preallocate
            SuperRat(ses).units(j).PFexist=zeros(1,4);
            SuperRat(ses).units(j).FiresDuringRun=zeros(1,4);
            SuperRat(ses).units(j).RunRates=repmat({[nan]},1,4);
            SuperRat(ses).units(j).FieldProps=FieldProps;
            
            % if no valid spikes, move on
            if isempty(allspikes)
                fprintf('\nSess %d Unit %d (%d tot) has no valid spikes \n', ses,j,length(spikes));
                continue;
            end
            % initialize the field props in case this cell doesnt ever fire
            if runblocks
                %{
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
                %}
                %  now one for all the blocks
            else % if run whole day as unitary session
                for tr=1:4 % for each trajectory
                    % this this trajectory (e.g. start and finish are
                    % correct
                    keepinds=fulltraj(:,4)==trajinds(tr,1) & fulltraj(:,5)==trajinds(tr,2);
                    % did the animal spend enough time in this trajectory?
                    if sum(keepinds)>30*10 % 15 seconds of data
                        % pull this traj, order epochs chronologically
                        temptraj=sortrows(fulltraj(keepinds,:),1);
                        % 50 bin span, but a std of only 1/4 second
                        % or like 3/4 of a second
                        tooslow=SmoothMat2(temptraj(:,7),[0 50],veltimesmooth)<=speedthreshold;
                        temptraj(tooslow,:)=[];
                        % find the epochs so you can kill bad spikes
                        breaks=find(diff(temptraj(:,1))>1); % only use epochs that last more than a seconds
                        runepochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1); temptraj(end,1)]];
                        
                        runepochs(diff(runepochs,1,2)<=.5,:)=[]; 
                        temptraj=EpochCoords(temptraj,runepochs); % remove the really short runs
                        
                        % grab all the spikes WITHIN the run epochs
                        [Espikes,Erates,~,~,trspikes]=event_spikes(spikes(:,1),runepochs(:,1),0,diff(runepochs,1,2));
                        SuperRat(ses).units(j).RunRates{1,tr}=Erates; % run the splitter score later
                        % filter spike and lfp adata below to match coords
                        
                        % we ahve the rates at each run, theyre zero, so
                        % now we cant build place maps
                        if length(Espikes)<20
                            %fprintf('no spikes detected ses %d unit %d tr %d (%.2f secs) \n',...
                            %ses, j, tr, sum(keepinds)/30);
                            continue;% if not enough spikes, go to next trajectory
                        end
                        
                        spikepos=interp1(temptraj(:,1),temptraj(:,8),Espikes,'nearest'); % now interp to position
                        % now you can calculate precession and phase
                        % preference
                        if calctheta
                            spikephase=interp1(eegmat(:,1),eegmat(:,3),Espikes,'nearest');
                            thetastats.spikeprefP(tr)=circ_rtest(spikephase);
                            thetastats.spikeprefV(tr)=circ_r(spikephase);
                            thetastats.spikeprefrho(tr)=circ_mean(spikephase);
                            [thetastats.slope(tr),thetastast.phasestart(tr),~,thetastast.precessionpP]=...
                                corrC2Lin_Kempter2012(spikepos,spikephase);
                        end
                        
                        % get occupancy map
                        [bincts,bins]=histcounts(temptraj(:,8),1:100); % in bins
                        occupancy=bincts/30; % convert to real time (seconds)
                        smoothoccup=SmoothMat2(occupancy,[5 0],2); % smooth over 2 pixels
                        SuperRat(ses).LinOccup(tr,:)=smoothoccup;
                        fullcurves={}; cellsort=[];
                        
                        if verbose, figure; end
                        % now for each unit capture a mean tuning curve for each run
                        % skaggs is run on presmoothed spike and occcupancy
                        if runboots>0
                            if~isempty(Espikes)
                                Bmax=nan(1,runboots); Binfo=nan(1,runboots); Bsparsity=nan(1,runboots);
                                for bt=1:runboots
                                    EBspikes=trspikes;
                                    for i=1:length(trspikes)
                                        % circ shift the spikes in this epoch by a min of 1/2 second
                                        if ~isempty(trspikes{i})
                                            EBspikes{i}=PermuteSpikeTimes(trspikes{i},runepochs(i,:)',0.2);
                                        end
                                    end
                                    Bspikepos=interp1(temptraj(:,1),temptraj(:,8),cell2mat(EBspikes'),'nearest');
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
                        
                        %%%%%%%%%%
                        %calculate all the fieldprops
                        %%%%%%%%%%
                        % now qualify this as a pf or not
                        [FieldProps.PFmax(1,tr),FieldProps.PFmaxpos(1,tr)]=max(LinPlaceField); % has to be above ??
                        FieldProps.Zpfmax(1,tr)=max(nanzscore(LinPlaceField)); % has to be above ??
                        
                        % run skaggs on individually smoothed occupancy and
                        % spikes
                        [FieldProps.info(1,tr),FieldProps.sparsity(1,tr)]=Skaggs_basic(smoothoccup./max(smoothoccup),...
                            smoothspikes,nanmean(smoothspikes));
                        
                        % and now pf size e.g. size of top 75%% of max rate
                        % this will underestimate, because it finds the
                        % first that drops below the high rate 'flood
                        % fill' method
                        pfstart=find(LinPlaceField(1:FieldProps.PFmaxpos(1,tr)) < FieldProps.PFmax(1,tr)*.25,1,'last');
                        if isempty(pfstart), pfstart=1; end
                        pfend=FieldProps.PFmaxpos(1,tr)+find(LinPlaceField(FieldProps.PFmaxpos(1,tr):end)...
                            < FieldProps.PFmax(1,tr)*.25,1,'first');
                        if isempty(pfend), pfend=100; end
                        
                        % try to calculate pf size, otherwise it doesnt exist
                        FieldProps.PFsize(1,tr)=pfend-pfstart;
                        
                        % now add to the session struct
                        SuperRat(ses).units(j).LinPlaceFields{1}(tr,:)=LinPlaceField;
                        
                        
                        if runboots>0
                            FieldProps.PFmaxP(1,tr)=1-normcdf(FieldProps.PFmax(1,tr),nanmean(Bmax),nanstd(Bmax)); % peak statistically high?
                            FieldProps.infoP(1,tr)=1-normcdf(FieldProps.info(1,tr),nanmean(Binfo),nanstd(Binfo)); % information statistically high?
                            FieldProps.sparsityP(1,tr)=1-normcdf(FieldProps.sparsity(1,tr),nanmean(Bsparsity),nanstd(Bsparsity)); % Sparsity?
                            SuperRat(ses).units(j).PFexist(1,tr)=FieldProps.PFmaxP(1,tr)<.05...
                                && FieldProps.PFmax(1,tr)>2 && FieldProps.PFsize(1,tr)<75; % is it a place cell
                        else
                            SuperRat(ses).units(j).PFexist(1,tr)=FieldProps.PFmax(1,tr)>2 && FieldProps.Zpfmax(1,tr)>2 && FieldProps.PFsize(1,tr)<75; % is it a place cell
                        end
                        SuperRat(ses).units(j).FiresDuringRun(1,tr)=FieldProps.PFmax(1,tr)>1; % is this cell active? e.g. does it pass 1 hz at anywhere on track
                        
                        
                        
                        if verbose % defunct because i folded pfmax etc into fieldprops
                            drawnow;
                            title(sprintf('Info: %.2f P=%.5f Peak: %.2f, p=%.5f',FieldProps.info(tr),...
                                FieldProps.infoP(tr),FieldProps.PFmax(tr),FieldProps.PFmaxP(tr)));
                            answer = questdlg('Stop verbose?', 'Verbose Menu', 'Yes, quiet','no, more figs','Yes, quiet');
                            verbose=contains(answer,'no');
                        end
                    else
                        fprintf('Sess %d traj %d doesnt have enough runs (%s sec)',ses,tr,sum(keepinds)/30);
                    end % if we have trajectories
                    
                end % of trajectories
            end % of block breakout
            if calctheta, SuperRat(ses).units(j).thetastats=thetastats; end
            SuperRat(ses).units(j).FieldProps=FieldProps; % add all the props
            fprintf('%d ',j);
        end % units
        fprintf('\n Session %d done in %.2f mins \n',ses,toc/60);
    end
end
if any(savedir~=0)
    save(fullfile(savedir,'ClaireData4-3-20'),'SuperRat');
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
        if length(SuperRat(i).units(j).RunRates)>=2
            if ~isempty(SuperRat(i).units(j).RunRates{1}) && ~isempty(SuperRat(i).units(j).RunRates{2})
                % bl takes care if youre breaking out by block
                
                RunCell=SuperRat(i).units(j).RunRates;
                RunRates=[RunCell{1} RunCell{2};...
                    ones(1,length(RunCell{1})) ones(1,length(RunCell{2}))+1]';
                % out left vs out right (1, 2), in left vs in right
                % positive is left preferring, negative is right preferring
                % for now these scores are total spikes, maybe well need max
                % which would take a different math here for the P value
                TrajScore=(nanmean(RunRates(RunRates(:,2)==1,1))-nanmean(RunRates(RunRates(:,2)==2,1)))/...
                    (nanmean(RunRates(RunRates(:,2)==1,1))+nanmean(RunRates(RunRates(:,2)==1,1)));
                
                if boot>0
                    % get the pull for
                    outTraj=RunRates(RunRates(:,2)==1 | RunRates(:,2)==2,1);
                    
                    TrajBoot=[nan nan];
                    for bt=1:boot
                        % randomize the order so you can even odd them
                        outOrder=randperm(length(outTraj));
                        TrajBoot(bt,1)=abs(nanmean(outTraj(mod(outOrder,2)==1))-...
                            nanmean(outTraj(mod(outOrder,2)==0)))/...
                            (nanmean(outTraj(mod(outOrder,2)==1))+...
                            nanmean(outTraj(mod(outOrder,2)==0)));
                    end
                    TrajScore(2)=1-normcdf(abs(TrajScore(1)), nanmean(TrajBoot(:,1)),nanstd(TrajBoot(:,1)));
                    
                else
                    TrajScore(2)=nan;
                end
                SuperRat(i).units(j).TrajScores=TrajScore;
            end
            waitbar(j/length(SuperRat(i).units),ha,sprintf('running unit %d',j));
        end
    end
    close(ha);
    fprintf('Session %d took %.2f minutes \n',i,toc/60);
end





%% just look at some cells to see if they're selective and that
% the raw rates look reasonable

ses=1;
% now lets see if there are valid priors here
for i=1:length(SuperRat(ses).units)
    if any(SuperRat(ses).units(i).PFexist)
        figure;
        try, plot(SuperRat(ses).units(i).LinPlaceFields{1}(1,:)); end
        hold on;
        try, plot(-SuperRat(ses).units(i).LinPlaceFields{1}(2,:)); end
        try, plot(SuperRat(ses).units(i).LinPlaceFields{1}(3,:)); end
       try,  plot(-SuperRat(ses).units(i).LinPlaceFields{1}(4,:)); end
        
        title([num2str(i) ' ' SuperRat(ses).units(i).area SuperRat(ses).units(i).type]);
        
    end
end
% looks really good actually!

%%



