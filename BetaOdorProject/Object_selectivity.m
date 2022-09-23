% now lets add some object selectivity onto this

% first, how many cells show object selectivity, lets see if we have the
% data in the behavior first


%load('E:\Claire Data\ClaireDataFull-LinPos-LongTrack.mat')


% object epochs are 'sniff' and the right and left 'correct' are the sides,
% but if its inco0rrect, you have to flip the odor

% cells that were active during the recording are separated first, these
% cells have to fire during the 'run' sessions, so basically they have to
% fire spikes during the trial period (take the grand sum of trial start
% and ends)

% odor selective is whether the cell fires more for one odor over the other
% (using a bootstrap with p<.05)

% taskresponsive is whether they exhibited significant changes during odor
% period (compared to a trial matched pre odor period) and this is a
% signrank test by trial

%% This plots the left and right rate for each cell locked to odor onset
% these are single cell example files

for ses=1:length(SuperRat)
    
    timeedges=-2:.001:2; timebins=-2:.001:1.999;
    % so maybe something like a left plot mean rate patch, then right plot
    % a colormap of abs diff/sum?
    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    mytrialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10];
    mytrialmat(mytrialmat(:,4)==0,:)=[]; %remove incorrect trials
    [temp,sortinds]=sortrows(mytrialmat,[3 1]); cutspot=find(diff(mytrialmat(sortinds,3)~=0)); odorid=mytrialmat(:,3);
    % okay now get the perievent raster, lets do the buzsaki way thats a
    % 100 msec smoothing kernel for each trial and grab -2 seconds to +2
    % seconds and patch the mean curve
    
    for i=1:length(SuperRat(ses).units)
        
        % get spikes
        [~,~,~,~,~,spikets]=event_spikes(SuperRat(ses).units(i).ts,mytrialmat(:,1),abs(timeedges(1)),timeedges(end));
        % now smooth each into a vector
        startcurves=cellfun(@(a) SmoothMat2(histcounts(a,timeedges),[1000 0],100), spikets, 'UniformOutput', false);
        
        % or unsmoothed
        %curves=cellfun(@(a)histcounts(a,-2:.001:2), spikets, 'UniformOutput', false);
        curvemat=cell2mat(startcurves');
        
        % now get the trial rates for object selectivity
        [~,trspikes]=event_spikes(SuperRat(ses).units(i).ts,mytrialmat(:,1),0,mytrialmat(:,2)-mytrialmat(:,1));
        % and now get the pre trial rates for a trial-length matched control
        [~,pretrspikes]=event_spikes(SuperRat(ses).units(i).ts,mytrialmat(:,1)-(mytrialmat(:,2)-mytrialmat(:,1)),0,mytrialmat(:,2)-mytrialmat(:,1));
        
        coderp=ranksum(trspikes(mytrialmat(:,3)==1),trspikes(mytrialmat(:,3)==0));
        responderp=signrank(trspikes,pretrspikes);
        if coderp<.05 || responderp<.05
            figure;
            sp(1)=subplot(2,1,1);
            imagesc(timebins, 1:size(startcurves,2), curvemat(sortinds,:)); hold on;
            plot([timebins(1) timebins(end)],[cutspot cutspot],'k','LineWidth',2);
            sp(2)=subplot(2,1,2);
            plot(timebins,nanmean(curvemat(odorid==1,:)),'r','LineWidth',3);
            
            patch([timebins fliplr(timebins)],[nanmean(curvemat(odorid==1,:))+nanstd(curvemat(odorid==1,:),0,1)...
                fliplr(nanmean(curvemat(odorid==1,:))-nanstd(curvemat(odorid==1,:),0,1))],[.85 .325 .098],...
                'LineStyle','none','FaceAlpha',.5)
            hold on; plot(timebins,nanmean(curvemat(odorid==0,:)),'b','LineWidth',3);
            patch([timebins fliplr(timebins)],[nanmean(curvemat(odorid==0,:))+nanstd(curvemat(odorid==0,:),0,1)...
                fliplr(nanmean(curvemat(odorid==0,:))-nanstd(curvemat(odorid==0,:),0,1))],[0    0.4470    0.7410],...
                'LineStyle','none','FaceAlpha',.5)
            linkaxes(sp,'x');
            % now get the overall rate
            title(sprintf('Mean Sampling period = %.2f P for responsive is %.2f p for selective is %.2f',...
                mean(mytrialmat(:,2)-mytrialmat(:,1)),responderp,coderp));
        end
    end
end

%%
% this calculates the following:
% first all active cells (those who fire at all during odor period)
% second all odor responsive cells (chg pre to post odor onset)
% third odor selective cells (odor 1 vs odor 2 only correct trials)
% fourth odor locked tuning curves
% fifth 

useblocks=false;
timeLock='end';
maxTime=2;
bootct=1000;
pCrit=1-normcdf(1.5); % noted in paper!!!!
timeedges=-2:.001:2; timebins=-2:.001:1.999;

% generate a repeatable random seed
rstream = RandStream('dsfmt19937','Seed',16);
RandStream.setGlobalStream(rstream);

% this adds the following fields:


% OdorRates: a nx3 matrix, 1 rate, 2 lr10, 3 ci10
% odorSelective: a 4 column table, 1 diff in mean rates, 2 pval, 3 is
%   p<05, 4 is the dprime effect size
%   second row is for incorrect trials
% taskResponsive: odor rate, preodor rate, pvalue
% Curves

for ses=1:length(SuperRat)

    tic
    wb=waitbar(0,'Starting to run cells');

    %%%%% grab trialmat
    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    % mat is 1 start, 2 end 3 left right 4 correct incorrect and 5 epoch
    fulltrialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
    % trials have to be in an analyzed block, and have to be >.5 seconds
    trialUse=ismember(fulltrialmat(:,5),SuperRat(ses).RunEpochs) & ...
        (fulltrialmat(:,2)-fulltrialmat(:,1))>=.5;

    trialmat=fulltrialmat(trialUse,:); % only take correct trials
    trialmat(:,6)=trialmat(:,2)-trialmat(:,1); % full poke time
    % if we dont want to use the full poke time
    if ~isnan(maxTime) 
        trialmat=trialmat(trialmat(:,6)<maxTime,:);
    end
    

    %%%%%% grab  fields
    % odorRates- table of rates by trial, rates, then odor (lr10), then
    % ci10
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorRates'); end
    % odorSelective- score (diff/sum), pval, is sig, dprime effect size
    % first row is correct, second is incorrect
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorSelective'); end
    % odorMeans- unnecessary
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorMeans'); end
    % odorResponsive- unnecessary
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorResponsive'); end
    % task responsive- table, odor rate, preodor rate, pval of signrank
    % test
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'taskResponsive'); end
    % curves- vector by time, selectivity index for correct then incorret
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'curves'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'startcurves','endcurves'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'activeCell'); end

    %%%%% for each cell...
    for i=1:length(SuperRat(ses).units)
        %%%%% initiate variables
         lrmeans=table('Size',[2 2],'VariableTypes',repmat({'double'},1,2),...
            'VariableNames',{'Left','Right'},'RowNames',{'Corr','Incorr'});
        Selectivitydata=table('Size',[2 4],'VariableTypes',repmat({'double'},1,4),...
            'VariableNames',{'score','pval','issig','dprime'},'RowNames',{'Corr','Incorr'});
        startcurves=[{},{}]; endcurves=[{},{}];
        taskResponsive=nan(1,3); % before, after, pval of signrank
        Selectivitydata.score(1)=nan;

        % first see if cell even fires
        transitions=diff(ismember([0; SuperRat(ses).tracking.data(1:end-1,6);0],SuperRat(ses).RunEpochs));
        epochTimes=[SuperRat(ses).tracking.data(transitions>0,1) SuperRat(ses).tracking.data(transitions<0,1)];
        [nspks,~,blockspikes]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
                epochTimes(:,1),0,diff(epochTimes,1,2));
        

        % gather the spkevs and trialspikes
         [allspks,allspkevs,~,trspks]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
                trialmat(:,1),0,trialmat(:,2)-trialmat(:,1));
         
         
         % remove blocks that arent used, or dont include enough spikes
         if useblocks
            % this removes blocks of trials where the cell doesnt fire,
            % this is somewhat of a bandaid if the cell isnt stable, so I
            % dont use this

            % remove trials from blocks where the cell has fewer spikes than
            % trials (I dont think she does this)
            spknums=cellfun(@(a) length(a), blockspikes);
            spikesperblock=accumarray(trialmat(:,5),spknums);
            trialsperblock=accumarray(trialmat(:,5),1);
            %find those inds and remove them
            keepblocktrials=sum(find(spikesperblock>trialsperblock/2)'==trialmat(:,5),2)>0;
            mytrialmat=trialmat(keepblocktrials,:);
        else
            mytrialmat=trialmat;

         end
       
         % if length(nspks)>100
        if length(allspks)>=10 %length(mytrialmat) % cell has to have more odor spikes than there were trials


            % grab odors, and grab correct
            odorid=mytrialmat(:,3); isCorr=mytrialmat(:,4)==1;
            % get the spike rate vectors
            [totalspikes,spkevs,~,spikeinds]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
                mytrialmat(:,1),0,mytrialmat(:,2)-mytrialmat(:,1)); % 0 to at least .5 secs
            [~,prespkevs]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
                mytrialmat(:,1),mytrialmat(:,2)-mytrialmat(:,1),0); % matching at least -.5 to 0
            
            % this is mean rate at odor, mean rate before odor, signrank p
            % value
            taskResponsive=[mean(spkevs(isCorr),'omitnan') mean(prespkevs(isCorr),'omitnan'),...
                signrank(spkevs(isCorr),prespkevs(isCorr))];
            
            % get the mean rates for each
            % LR10 left first, right second
            for k=1:2 % 1 is correct, 2 is incorrect
                lrmeans.Left(k)=mean(spkevs(odorid==1 & isCorr==2-k),'omitnan'); % and this is spikes per second
                lrmeans.Right(k)=mean(spkevs(odorid==0 & isCorr==2-k),'omitnan');
                %%%%% selectivitydata

                Selectivitydata.score(k)=diff(lrmeans{k,:})/sum(lrmeans{k,:}); % 
                % and selectivity index
                [~,Selectivitydata.pval(k)]=SelectivityIndex(spkevs(isCorr==2-k),odorid(isCorr==2-k),bootct); % 500 boots
                Selectivitydata.issig(k)=Selectivitydata.pval(k) <= pCrit & ~isnan(Selectivitydata.score(k)); % boolean
                Selectivitydata.dprime(k)=dprime(spkevs(odorid==1 & isCorr==2-k),spkevs(odorid==0 & isCorr==2-k)); % and a dprime

                % Now create the curves
                 [~,~,~,~,~,spikets]=event_spikes(SuperRat(ses).units(i).ts,...
                        trialmat(:,1),abs(timeedges(1)),timeedges(end));
             
                % now smooth each into a rate vector for correct and incorrect
                % I only need to correct for ntrials (cause bins are the
                % same
                rawCurves=SmoothMat2(histcounts(cell2mat(spikets(odorid==1 & isCorr==2-k)'),...
                    timeedges)./sum(odorid==1 & isCorr==2-k),[1000 0],50);
                rawCurves(2,:)=SmoothMat2(histcounts(cell2mat(spikets(odorid==0 & isCorr==2-k)'),...
                    timeedges)./sum(odorid==0 & isCorr==2-k),[1000 0],50);
                % this is the SI value here for k=1 correct and k=2
                % incorrect
                startcurves{k}=diff(rawCurves)./sum(rawCurves);
                % repeat for withdrawal locked curves (always in msec,
                % always centered around 0
                 [~,~,~,~,~,spikets]=event_spikes(SuperRat(ses).units(i).ts,...
                        trialmat(:,2),abs(timeedges(1)),timeedges(end));
                rawCurves=SmoothMat2(histcounts(cell2mat(spikets(odorid==1 & isCorr==2-k)'),timeedges),[1000 0],50);
                rawCurves(2,:)=SmoothMat2(histcounts(cell2mat(spikets(odorid==0 & isCorr==2-k)'),timeedges),[1000 0],50);
                endcurves{k}=diff(rawCurves)./sum(rawCurves);

            end
            SuperRat(ses).units(i).activeCell=true;
        else
            SuperRat(ses).units(i).activeCell=false;
        end
        SuperRat(ses).units(i).taskResponsive=taskResponsive;
        %SuperRat(ses).units(i).OdorResponsive=odorResponse; % this is not
        %really useful, taskResponsive is better
        SuperRat(ses).units(i).OdorRates=[allspkevs, trialmat(:,3:4)];
        SuperRat(ses).units(i).OdorMeans=lrmeans; % upload to the struct
        SuperRat(ses).units(i).OdorSelective=Selectivitydata;
        SuperRat(ses).units(i).startcurves=startcurves;
        SuperRat(ses).units(i).endcurves=endcurves;
        
        waitbar(i/length(SuperRat(ses).units),wb,sprintf('Running unit %d',i));
    end
    close(wb);
    % record how long it took (it usually is super short)
    fprintf('Ses %d took %.2f seconds \n',ses,toc);
    
end



%% Repeat for reward poke firing

%************ DONT HAVE ALL THE REWARD DATA THOUGH ****************

%{
% this is claires way
useblocks=false;
bootct=500;

for ses=1:length(SuperRat)
    tic
    % so maybe something like a left plot mean rate patch, then right plot
    % a colormap of abs diff/sum?
    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    % mat is 1 start, 2 end 3 left right 4 correct incorrect and 5 epoch\
    if isfield(trialdata,'rewardstart')
        mytrialmat=[trialdata.rewardstart trialdata.rewardend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
        [temp,sortinds]=sortrows(mytrialmat,[3 1]); odorid=mytrialmat(:,3);
        % okay now get the perievent raster, lets do the buzsaki way thats a
        % 100 msec smoothing kernel for each trial and grab -2 seconds to +2
        % seconds and patch the mean curve
        rlmeans=nan(length(SuperRat(ses).units),2);
        wb=waitbar(0,'Starting to run cells');
        try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'RewardRates'); end
        try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'RewardMeans'); end
        try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'RewardSelective'); end
        
        for i=1:length(SuperRat(ses).units)
            
            %build  the spikes
            [~,spkevs,~,~,~,controls]=event_spikes(SuperRat(ses).units(i).ts(:,1),mytrialmat(:,1),0,mytrialmat(:,2)-mytrialmat(:,1));
            %curves=cellfun(@(a) SmoothMat2(histcounts(a,0:.001:1),[1000 0],100), controls, 'UniformOutput', false);
            %curvemat=cell2mat(curves');
            if useblocks
                epochcounts=unique(mytrialmat(:,5)); % how many sessions are in this recording day?
                Selectivitydata=[]; rlmeans=[];
                for ep=1:length(unique(epochcounts)) % calc for each session
                    % get all the spikes for the correct trials and this epoch
                    goodtrials=mytrialmat(:,4)==1 & mytrialmat(:,5)==epochcounts(ep); %
                    if nanmean(any(spkevs(goodtrials)))>=0.5 % if the cell spikes on at least half the trials
                        % get the mean rates for each
                        % LR10 left first, right second
                        rlmeans(ep,1)=nanmean(spkevs(goodtrials & mytrialmat(:,3)==1)); % and this is spikes per second
                        rlmeans(ep,2)=nanmean(spkevs(goodtrials & mytrialmat(:,3)==0));
                        % and an effect size just to see
                        %rlmeans(ep,3)=dprime(spkevs(usetrials & trialmat(:,3)==1),spkevs(usetrials & trialmat(:,3)==0));
                        % and selectivity index
                        [Selectivitydata(ep,1),Selectivitydata(ep,2)]=SelectivityIndex(spkevs(goodtrials),mytrialmat(goodtrials,3),bootct); % 500 boots
                        Selectivitydata(ep,3)=Selectivitydata(ep,2) <= 0.05 & ~isnan(Selectivitydata(ep,1)); % how many boots
                    else
                        rlmeans(ep,1:3)=nan; Selectivitydata(ep,1:3)=nan; % if not enouch spikes, its nanned out
                    end
                    SuperRat(ses).units(i).RewardRates{ep}=[spkevs(goodtrials) mytrialmat(goodtrials,3)];
                    SuperRat(ses).units(i).RewardMeans(ep,:)=rlmeans; % upload to the struct
                    SuperRat(ses).units(i).RewardSelective(ep,:)=Selectivitydata;
                end
            else
                goodtrials=mytrialmat(:,4)==1;  % only correct trials
                if nanmean(any(spkevs))>=0.5 % if the cell spikes on at least half the trials
                    % get the mean rates for each
                    % LR10 left first, right second
                    rlmeans=nanmean(spkevs(goodtrials & mytrialmat(:,3)==1)); % and this is spikes per second
                    rlmeans(1,2)=nanmean(spkevs(goodtrials & mytrialmat(:,3)==0));
                    % and an effect size just to see
                    %rlmeans(ep,3)=dprime(spkevs(usetrials & trialmat(:,3)==1),spkevs(usetrials & trialmat(:,3)==0));
                    % and selectivity index
                    [Selectivitydata,Selectivitydata(1,2)]=SelectivityIndex(spkevs(goodtrials),mytrialmat(goodtrials,3),bootct); % 500 boots
                    Selectivitydata(1,3)=Selectivitydata(2) <= 0.05 & ~isnan(Selectivitydata(1)); % how many boots
                else
                    rlmeans(1,1:3)=nan; Selectivitydata(1,1:3)=nan; % if not enouch spikes, its nanned out
                end
                SuperRat(ses).units(i).RewardRates=[spkevs(goodtrials)' mytrialmat((goodtrials),3)];
                SuperRat(ses).units(i).RewardMeans=rlmeans; % upload to the struct
                SuperRat(ses).units(i).RewardSelective=Selectivitydata;
            end
            
            waitbar(i/length(SuperRat(ses).units),wb,sprintf('Running unit %d',i));
        end
        close(wb);
        % record how long it took (it usually is super short)
        fprintf('Ses %d took %.2f seconds \n',ses,toc);
    else
        fprintf('Ses %d doesnt have reward events \n',ses);
        SuperRat(ses).units(:).RewardRates=nan(1,3);
        
    end
    
end
%}

%% there is a discrepancy basically and i dont remember how we reconsiled it....
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% analysis starts here




%% here is my total number of pyrams and INS that are responsive and selective
varTypes=repmat({'double'},1,7);
responseTable=table('size',[2 7],'VariableTypes',varTypes,'VariableNames',{'allTot','pyrTot','pyrResp','pyrSel','inTot','inResp','inSel'},...
    'RowNames',{'PFC','CA1'});
allcells=cell2mat({SuperRat.units});
allcells=allcells([allcells.activeCell]==1);
regions={'PFC','CA1'};
for i=1:2
    Pcells=allcells(strcmpi({allcells.area},regions{i}));
    responseTable.allTot(i)=length(Pcells);
    Ppyrams=Pcells(strcmpi({Pcells.type},'pyr'));
    responseTable.pyrTot(i)=length(Ppyrams); % total number of cells
    responseTable.pyrResp(i)=sum(cellfun(@(a) a(3)<.05, {Ppyrams.taskResponsive})); % of tot that are task responsive

    % now % of task responsive that are odor selective
    responseTable.pyrSel(i)=sum(cellfun(@(a) a{1,3}==1, {Ppyrams(cellfun(@(a) a(3)<.05, {Ppyrams.taskResponsive})).OdorSelective}));

    Pins=Pcells(cellfun(@(a) contains(a,'in'), {Pcells.type}));
    responseTable.inTot(i)=length(Pins);
    responseTable.inResp(i)=sum(cellfun(@(a) a(3)<.05, {Pins.taskResponsive}));
    responseTable.inSel(i)=sum(cellfun(@(a) a{1,3}==1, {Pins(cellfun(@(a) a(3)<.05, {Pins.taskResponsive})).OdorSelective}));
end

pctTable=responseTable;
for i=1:2
    pctTable{i,2}=responseTable{i,2}/responseTable{i,1};
    pctTable{i,3}=responseTable{i,3}/responseTable{i,2};
    pctTable{i,4}=responseTable{i,4}/responseTable{i,3};
    pctTable{i,5}=responseTable{i,5}/responseTable{i,1};
    pctTable{i,6}=responseTable{i,6}/responseTable{i,5};
    pctTable{i,7}=responseTable{i,7}/responseTable{i,6};
end



openvar('responseTable');
openvar('pctTable')

%% now a retabulation of all the task responsive and task selective units
% claire grossly overestimated the number of task responsive and task
% selective units, and its unclear why... one possibility is that she used
% a strikingly low threshold of 1.5 sd above the null distribution




%% I need to redo claires figure 2e and i- will need to make Selective a table
% also get the rates asa table
% need to recalc f and J and also show that SI for incorrect is weaker than
% for correct

% okay so we can reshape figure 2a and b using the above tables
% responsetable and pctTable

pcrit=1-normcdf(1.5);
pcrit=.05;

% and to get e and I,...
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
type={'pyr','in'};
load redToBlue; % the red to blue colormap
for t=1:2
    for i=1:2
        regcoders=allcells(strcmpi({allcells.area},regions{i}) &...
            cellfun(@(a) a{1,2}<pcrit,{allcells.OdorSelective}) &...
            strcmpi({allcells.type},type{t}));
        corrCurves=cell2mat(cellfun(@(a) a{1}, {regcoders.startcurves},'UniformOutput',false)');
        incCurves=cell2mat(cellfun(@(a) a{2}, {regcoders.startcurves},'UniformOutput',false)');
        % 2000 to 3000 is 0 to 1 second following the flag start
        [~,rowsort]=sort(mean(-corrCurves(:,2000:3000),2));
        figure('Position',[300+200*i+50*t 300-50*t 600 500]);

        subplot(2,3,1);
        imagesc(0:.001:1,1:size(corrCurves,1),corrCurves(rowsort,2000:3000))
        set(gca,'Colormap',redToBlue);
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);

        subplot(2,3,4);
        imagesc(0:.001:1,1:size(incCurves,1),incCurves(rowsort,2000:3000))
        set(gca,'Colormap',redToBlue); xlabel(sprintf('seconds from \n odor start'));
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);

        % and locked to exit
        corrCurves=cell2mat(cellfun(@(a) a{1}, {regcoders.endcurves},'UniformOutput',false)');
        incCurves=cell2mat(cellfun(@(a) a{2}, {regcoders.endcurves},'UniformOutput',false)');
        % 2000 to 3000 is 0 to 1 second following the flag start
        [~,rowsort]=sort(mean(-corrCurves(:,1000:2000),2));

        subplot(2,3,2);
        imagesc(-1:.001:0,1:size(corrCurves,1),corrCurves(rowsort,1000:2000))
        set(gca,'Colormap',redToBlue);
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);

        subplot(2,3,5);
        imagesc(-1:.001:0,1:size(incCurves,1),incCurves(rowsort,1000:2000))
        set(gca,'Colormap',redToBlue); xlabel(sprintf('seconds from \n odor exit'))
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);


        SIvals=cell2mat(cellfun(@(a) a.score, {regcoders.OdorSelective},'UniformOutput',false));
        %Weights=cell2mat(cellfun(@(a) -log10(a.pval(1)), {regcoders.OdorSelective},'UniformOutput',false));
        %Weights(Weights==inf | Weights==-inf)=5; % this happens if its beyond precision i think...
        sp=subplot(3,3,6);
        %mdl=fitlm(SIvals(1,:),SIvals(2,:),'Weights',Weights);
        mdl=fitlm(SIvals(1,:),SIvals(2,:));
        [r2]=corr(SIvals(1,:)',SIvals(2,:)','rows','complete');
        plot(sp,mdl); kids=get(gca,'Children');
        set(kids(4),'Marker','none'); hold on;
        %scatter(SIvals(1,:),SIvals(2,:),rescale(Weights,10, 32),colors(i,:),'filled');
        scatter(SIvals(1,:),SIvals(2,:),16,colors(i,:),'filled');
        kids=get(gca,'Children'); legend(kids([1 3 4]))
        title(sprintf('Slope %.2f \n P %.2e',r2,mdl.Coefficients.pValue(2)));
        legend off;
        for cl=2:4, set(kids(cl),'color','k'); end
        %scatter(SIvals(1,:),SIvals(2,:),10,colors(i,:),'filled');
        axis tight; xlabel('Correct SI'); ylabel('Incorrect SI');
         box off;
        sgtitle([regions{i} ' ' type{t}]);
    end
end

%% and for noncoder pyrams
for t=1:2
    for i=1:2
        regcoders=allcells(strcmpi({allcells.area},regions{i}) &...
            cellfun(@(a) a{1,2}>pcrit,{allcells.OdorSelective}) &...
            cellfun(@(a) a(3)<.05,{allcells.taskResponsive}) & ...
            strcmpi({allcells.type},type{t}));
        corrCurves=cell2mat(cellfun(@(a) a{1}, {regcoders.startcurves},'UniformOutput',false)');
        incCurves=cell2mat(cellfun(@(a) a{2}, {regcoders.startcurves},'UniformOutput',false)');
        % 2000 to 3000 is 0 to 1 second following the flag start
        [~,rowsort]=sort(mean(-corrCurves(:,2000:3000),2));
        figure('Position',[300+200*i+50*t 300-50*t 600 500]);

        subplot(2,3,1);
        imagesc(0:.001:1,1:size(corrCurves,1),corrCurves(rowsort,2000:3000))
        set(gca,'Colormap',redToBlue);
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);

        subplot(2,3,4);
        imagesc(0:.001:1,1:size(incCurves,1),incCurves(rowsort,2000:3000))
        set(gca,'Colormap',redToBlue); xlabel(sprintf('seconds from \n odor start'));
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);

        % and locked to exit
        corrCurves=cell2mat(cellfun(@(a) a{1}, {regcoders.endcurves},'UniformOutput',false)');
        incCurves=cell2mat(cellfun(@(a) a{2}, {regcoders.endcurves},'UniformOutput',false)');
        % 1000 to 2000 is 0 to 1 second before the flag
        [~,rowsort]=sort(mean(-corrCurves(:,1000:2000),2));

        subplot(2,3,2);
        imagesc(-1:.001:0,1:size(corrCurves,1),corrCurves(rowsort,1000:2000))
        set(gca,'Colormap',redToBlue);
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);

        subplot(2,3,5);
        imagesc(-1:.001:0,1:size(incCurves,1),incCurves(rowsort,1000:2000))
        set(gca,'Colormap',redToBlue); xlabel(sprintf('seconds from \n odor exit'))
        cmax=max(abs(get(gca,'CLim'))); set(gca,'CLim',[-cmax cmax]);


        SIvals=cell2mat(cellfun(@(a) a.score, {regcoders.OdorSelective},'UniformOutput',false));
        %Weights=cell2mat(cellfun(@(a) -log10(a.pval(1)), {regcoders.OdorSelective},'UniformOutput',false));
        %Weights(Weights==inf | Weights==-inf)=5; % this happens if its beyond precision i think...
        sp=subplot(3,3,6);
        %mdl=fitlm(SIvals(1,:),SIvals(2,:),'Weights',Weights);
        mdl=fitlm(SIvals(1,:),SIvals(2,:));
        [r2]=corr(SIvals(1,:)',SIvals(2,:)','rows','complete');
        plot(sp,mdl); kids=get(gca,'Children');
        set(kids(4),'Marker','none'); hold on;
        %scatter(SIvals(1,:),SIvals(2,:),rescale(Weights,10, 32),colors(i,:),'filled');
        scatter(SIvals(1,:),SIvals(2,:),16,colors(i,:),'filled');
        kids=get(gca,'Children'); legend(kids([1 3 4]))
        title(sprintf('Slope %.2f \n P %.2e P %.2e',r2,mdl.Coefficients.pValue(2)));
        legend off;
        for cl=2:4, set(kids(cl),'color','k'); end
        %scatter(SIvals(1,:),SIvals(2,:),10,colors(i,:),'filled');
        axis tight; xlabel('Correct SI'); ylabel('Incorrect SI');
        sgtitle(['Noncoders ' regions{i} ' ' type{t}]);
    end
end

% massive variable to clear
clear regcoders;
%% and now some dicing of what the hell is goign on with 'taskresponsive'

% this is a problem, a large proportion of 'task responsive' cells dont
% fire at all during the odor period.  I will have to cut these cells out
% back in the normal analysis.  This ought to only mess with a few
% statistics though....
for tp=1:2
    for i=1:2
        figure;
        % gather all task responsive cells
        responders=allcells(cellfun(@(a) a(3)<.05,{allcells.taskResponsive}) &...
            strcmpi({allcells.type},type{tp}) & strcmpi({allcells.area},regions{i}) &...
            [allcells.activeCell]);
        % tabulate the pre, post, and the diff
        % taskresponsive- 1. odor rate, 2. preodor rate, 3. signrank diff pval
        allrates=cell2mat({responders.taskResponsive}');
        
        scatter(allrates(:,1),allrates(:,1)-allrates(:,2));
        xlabel('odor rate'); ylabel('change (odor-preodor)')
        figure; histogram((allrates(:,1)-allrates(:,2))./(allrates(:,2)+allrates(:,1)),15);
        ylabel('counts'); xlabel('-1 preodor preferring <-> odor preferring 1');
        title(sprintf('%s, %s, %.1f percent drop rate', regions{i},type{tp},...
            mean(allrates(:,1)-allrates(:,2)<0)))

        figure;
        subplot(1,2,1);
        allBetaMVL=cell2mat({responders.betaMVL}');
        selectivities=cellfun(@(a) a.dprime(1), {responders.OdorSelective});
        scatter(allrates(:,1)-allrates(:,2),max(allBetaMVL,[],2));
        xlabel('chg in rate'); ylabel('betaMVL');
        subplot(1,2,2);
        scatter(allrates(:,1)-allrates(:,2),abs(selectivities));
        xlabel('chg in rate'); ylabel('abs(selectivity)');
    end
end
%%
% This is a remake of Claires Table S1
subcats={'Total','Active','Pyr','IN','TaskResp'};
myvars=['nSessions', strcat({'CA1 '},subcats), strcat({'PFC '},subcats)];
myrows={'1','2','3','4','5','6','7','8','Total'};
TableS1=table('Size',[9 11],'VariableTypes',repmat({'cell'},1,11),...
    'VariableNames',myvars,'RowNames',myrows);

% backwards of what i used to do
regions2={'CA1','PFC'};
[~,~,ratIDs]=unique({SuperRat.name});
TableS1.nSessions=num2str([accumarray(ratIDs,1); mean(accumarray(ratIDs,1))]);


for i=1:length(regions2)
    % first we do total:
    nTot=cellfun(@(a) sum(contains({a.area},regions2{i})), {SuperRat.units});
    ratTots=accumarray(ratIDs,nTot',[],@mean);
    nActive=cellfun(@(a) sum(contains({a.area},regions2{i}) &...
        [a.activeCell]), {SuperRat.units});
    ratActive=accumarray(ratIDs,nActive',[],@mean);
    nPyr=cellfun(@(a) sum(contains({a.area},regions2{i}) &...
        [a.activeCell] & contains({a.type},'pyr')), {SuperRat.units});
    ratPyr=accumarray(ratIDs,nPyr',[],@mean);   
    nInt=cellfun(@(a) sum(contains({a.area},regions2{i}) &...
        [a.activeCell] & contains({a.type},'in')), {SuperRat.units});
    ratInt=accumarray(ratIDs,nInt',[],@mean);
    nTaskresp=cellfun(@(a) sum(contains({a.area},regions2{i}) &...
        cellfun(@(b) b(3)<.05, {a.taskResponsive})), {SuperRat.units});
    ratTaskresp=accumarray(ratIDs,nTaskresp',[],@mean);


    for j=1:height(TableS1)-1
        TableS1.([regions2{i} ' ' subcats{1}])(j)={sprintf('%.2f',ratTots(j))};
        TableS1.([regions2{i} ' ' subcats{2}])(j)={sprintf('%.2f',ratActive(j))};
        TableS1.([regions2{i} ' ' subcats{3}])(j)={sprintf('%.2f',ratPyr(j))};
        TableS1.([regions2{i} ' ' subcats{4}])(j)={sprintf('%.2f',ratInt(j))};
        TableS1.([regions2{i} ' ' subcats{5}])(j)={sprintf('%.2f',ratTaskresp(j))};
    end
    % and totals
    TableS1.([regions2{i} ' ' subcats{1}])(end)={sprintf('%.2f +/-(%.2f)',...
        mean(ratTots), SEM(ratTots))};
    TableS1.([regions2{i} ' ' subcats{2}])(end)={sprintf('%.2f +/-(%.2f)',...
        mean(ratActive), SEM(ratActive))};
    TableS1.([regions2{i} ' ' subcats{3}])(end)={sprintf('%.2f +/-(%.2f)',...
        mean(ratPyr), SEM(ratPyr))};
    TableS1.([regions2{i} ' ' subcats{4}])(end)={sprintf('%.2f +/-(%.2f)',...
        mean(ratInt), SEM(ratInt))};
    TableS1.([regions2{i} ' ' subcats{5}])(end)={sprintf('%.2f +/-(%.2f)',...
        mean(ratTaskresp), SEM(ratTaskresp))};
end
openvar('TableS1')
%%
%
%
%
%
%
%    Reconsiliation work with Claires Accounting
%
%
%

%% lets get a list of each cells session and name and crossreference
% with claires tabulations


% reconsiliation work with Claire
%{
SuperUnits=orderfields(SuperRat(1).units);
%SuperUnits=rmfield(SuperUnits,'csi');
thisdate=repmat({[SuperRat(1).name ' day ' num2str(SuperRat(1).daynum)]},length(SuperUnits),1);

[SuperUnits.sess]=thisdate{:};
for i=2:length(SuperRat)
    theseunits=orderfields(SuperRat(i).units);
    % for some reason csi is in some of these structs
    if isfield(theseunits,'csi'), theseunits=rmfield(theseunits,'csi'); end
    thisdate=repmat({[SuperRat(i).name ' day ' num2str(SuperRat(i).daynum)]},length(theseunits),1);
    [theseunits.sess]=thisdate{:};
    SuperUnits=[SuperUnits theseunits];
end

% and the comparison to claires dataset
clairedata=load('E:\Brandeis datasets\Claire Data\ClaireObjResults.mat');
clairedata.names={'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
% so i have to remove all sessions with rat numbers
%  6 and 7

% okay now need to track down all the cell she thinks code...
clear ClaireRespCA1;
% first is her cell name, second is mine, 3 is p and 4 is is sig
for i=1:length(clairedata.npCellsCA1)
    ClaireRespCA1{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.npCellsCA1(i,1)},...
        clairedata.npCellsCA1(i,2),clairedata.npCellsCA1(i,3),...
        clairedata.npCellsCA1(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.npCellsCA1(i,1)},clairedata.npCellsCA1(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.npCellsCA1(i,3) & [SuperUnits.unitnum]==clairedata.npCellsCA1(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireRespCA1{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).taskResponsive;
        ClaireRespCA1{i,3}=odordata(2);
        ClaireRespCA1{i,4}=odordata(3);
        
    elseif clairedata.npCellsCA1(i,1)==6 || clairedata.npCellsCA1(i,1)==7
        ClaireRespCA1{i,2}='didnt analyze session';
    else
        ClaireRespCA1{i,2}='couldnt find cell';
    end
end


%okay now need to track down all the cell she thinks code...
clear ClaireSelCA1;
for i=1:length(clairedata.selectiveCA1)
    ClaireSelCA1{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.selectiveCA1(i,1)},...
        clairedata.selectiveCA1(i,2),clairedata.selectiveCA1(i,3),...
        clairedata.selectiveCA1(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.selectiveCA1(i,1)},clairedata.selectiveCA1(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.selectiveCA1(i,3) & [SuperUnits.unitnum]==clairedata.selectiveCA1(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireSelCA1{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).OdorSelective;
        ClaireSelCA1{i,3}=odordata(2);
        ClaireSelCA1{i,4}=odordata(3);
        
    elseif clairedata.selectiveCA1(i,1)==6 || clairedata.selectiveCA1(i,1)==7
        ClaireSelCA1{i,2}='didnt analyze session';
    else
        ClaireSelCA1{i,2}='couldnt find cell';
    end
end


% okay now need to track down all the cell she thinks code...
clear ClaireRespPFC;
for i=1:length(clairedata.npCellsPFC)
    ClaireRespPFC{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.npCellsPFC(i,1)},...
        clairedata.npCellsPFC(i,2),clairedata.npCellsPFC(i,3),...
        clairedata.npCellsPFC(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.npCellsPFC(i,1)},clairedata.npCellsPFC(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.npCellsPFC(i,3) & [SuperUnits.unitnum]==clairedata.npCellsPFC(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireRespPFC{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).taskResponsive;
        ClaireRespPFC{i,3}=odordata(2);
        ClaireRespPFC{i,4}=odordata(3);
        
    elseif clairedata.npCellsPFC(i,1)==6 || clairedata.npCellsPFC(i,1)==7
        ClaireRespPFC{i,2}='didnt analyze session';
    else
        ClaireRespPFC{i,2}='couldnt find cell';
    end
end


% okay now need to track down all the cell she thinks code...
clear ClaireSelPFC;
for i=1:length(clairedata.selectivePFC)
    ClaireSelPFC{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.selectivePFC(i,1)},...
        clairedata.selectivePFC(i,2),clairedata.selectivePFC(i,3),...
        clairedata.selectivePFC(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.selectivePFC(i,1)},clairedata.selectivePFC(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.selectivePFC(i,3) & [SuperUnits.unitnum]==clairedata.selectivePFC(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireSelPFC{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).OdorSelective;
        ClaireSelPFC{i,3}=odordata(2);
        ClaireSelPFC{i,4}=odordata(3);
        
    elseif clairedata.selectivePFC(i,1)==6 || clairedata.selectivePFC(i,1)==7
        ClaireSelPFC{i,2}='didnt analyze session';
    else
        ClaireSelPFC{i,2}='couldnt find cell';
    end
end

%}