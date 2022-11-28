 
%%
%{
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
%}

%% first do activeCells and taskResponsive

for ses=1:length(SuperRat)
    tic
    % get trial data
    trialdata=SuperRat(ses).trialdata;
    fulltrialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
    % trials have to be in an analyzed block, and have to be >.5 seconds
    trialUse=ismember(fulltrialmat(:,5),SuperRat(ses).RunEpochs) & ...
        (fulltrialmat(:,2)-fulltrialmat(:,1))>=.5 &...
        (fulltrialmat(:,2)-fulltrialmat(:,1))<2;
    trialmat=fulltrialmat(trialUse,:); % only take correct trials from high performing blocks

    odorid=trialmat(:,3); isCorr=trialmat(:,4)==1;


    % now get epochs
    transitions=diff(ismember([0; SuperRat(ses).tracking.data(1:end-1,6);0],SuperRat(ses).RunEpochs));
    epochTimes=[SuperRat(ses).tracking.data(transitions>0,1) SuperRat(ses).tracking.data(transitions<0,1)];
    

    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'taskResponsive'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'activeCell'); end

    for i=1:length(SuperRat(ses).units)
        % preallocate variables
        Selectivitydata=table('Size',[2 6],'VariableTypes',repmat({'double'},1,6),...
            'VariableNames',{'score','pval','issig','dprime','leftmean','rightmean'},'RowNames',{'Corr','Incorr'});
        startcurves=[{},{}]; endcurves=[{},{}];
        taskResponsive=nan(1,3); % before, after, pval of signrank
        Selectivitydata.score(1)=nan;

        % all spikes during each block
        [nspks,~,~,~,blockspikes]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
        epochTimes(:,1),0,diff(epochTimes,1,2));

        % gather the spkevs and trialspikes
        [allspks,spkevs,~,trspks]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
            trialmat(:,1),0,trialmat(:,2)-trialmat(:,1));

        % gather prespkevs for taskResponsive
        [~,prespkevs]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
            trialmat(:,1),trialmat(:,2)-trialmat(:,1),0); % matching at least -.5 to 0

        % all spikes in a +/- 2 second window surrounding odor sampling
        [allspks2]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
                trialmat(:,1),2,2);

        taskResponsive=nan(1,3); % before, after, pval of signrank

        %if length(nspks)>100 % active Units, 100 spikes across all run epochs
        %if length(nspks)>=length(mytrialmat) % cell has to have more odor spikes than there were trials
        %if size(SuperRat(ses).units(i).ts,1)>length(mytrialmat)
        if length(allspks2)>=length(trialmat)
        %if length(allspks)>=100


            % empiricallyif you 'runboth' and test on each odor
            % individually you actually get fewer cells that are task
            % responsive
            runboth=1;
            if runboth==1
                for o=1:2
                    taskResponsive(o,:)=[mean(spkevs(isCorr & odorid==o-1),'omitnan') mean(prespkevs(isCorr & odorid==o-1),'omitnan'),...
                        signrank(spkevs(isCorr & odorid==o-1),prespkevs(isCorr & odorid==o-1))];
                end
                [~,winner]=min(taskResponsive(:,3));
                taskResponsive=taskResponsive(winner,:);
            else
                taskResponsive=[mean(spkevs(isCorr),'omitnan') mean(prespkevs(isCorr),'omitnan'),...
                    signrank(spkevs(isCorr),prespkevs(isCorr))];
            end
            % first try left trials
            % then try right trials
            SuperRat(ses).units(i).activeCell=true;
        else
            SuperRat(ses).units(i).activeCell=false;
        end
        SuperRat(ses).units(i).taskResponsive=taskResponsive;
    end
    fprintf('session %d done, took %.2f seconds \n', ses,toc)
end

%% now odor rates...


% generate a repeatable random seed
rstream = RandStream('dsfmt19937','Seed',16);
RandStream.setGlobalStream(rstream);
bootct=1000;
odorPcrit=1-normcdf(1.5); % noted in paper!!!!
timeedges=-2:.001:2; timebins=-2:.001:1.999;
testing={'correct','incorrect'};

for ses=1:length(SuperRat)
    tic
    % get trial data
    trialdata=SuperRat(ses).trialdata;
    fulltrialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
    % trials have to be in an analyzed block, and have to be >.5 seconds
    trialUse=ismember(fulltrialmat(:,5),SuperRat(ses).RunEpochs) & ...
        (fulltrialmat(:,2)-fulltrialmat(:,1))>=.5 &...
        (fulltrialmat(:,2)-fulltrialmat(:,1))<2;
    trialmat=fulltrialmat(trialUse,:); % only take correct trials from high performing blocks
    odorid=trialmat(:,3); isCorr=trialmat(:,4)==1;

    % now get epochs
    transitions=diff(ismember([0; SuperRat(ses).tracking.data(1:end-1,6);0],SuperRat(ses).RunEpochs));
    epochTimes=[SuperRat(ses).tracking.data(transitions>0,1) SuperRat(ses).tracking.data(transitions<0,1)];


    mytrialmat=trialmat;
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorRates'); end
    % odorSelective- score (diff/sum), pval, is sig, dprime effect size
    % first row is correct, second is incorrect
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorSelective'); end
    % odorMeans- unnecessary
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorMeans'); end
    % odorResponsive- unnecessary
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorResponsive'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,{'startcurves','endcurves'}); end

    for i=1:length(SuperRat(ses).units)

        % preallocate variables
        Selectivitydata=array2table(nan(2,6),'VariableNames',{'score','pval','issig','dprime','leftmean','rightmean'},...
            'RowNames',{'Corr','Incorr'});

        startcurves=[{},{}]; endcurves=[{},{}];
        taskResponsive=nan(1,3); % before, after, pval of signrank


        % gather the spkevs and trialspikes
        [allspks,spkevs,~,trspks]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
            trialmat(:,1),0,trialmat(:,2)-trialmat(:,1));

        % spikes in trial matched pre odor period
        [~,prespkevs]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
            trialmat(:,1),trialmat(:,2)-trialmat(:,1),0); % matching at least -.5 to 0


        for k=1:2 % 1 is correct, 2 is incorrect
            try
                Selectivitydata.leftmean(k)=mean(spkevs(odorid==1 & isCorr==2-k),'omitnan'); % and this is spikes per second
                Selectivitydata.rightmean(k)=mean(spkevs(odorid==0 & isCorr==2-k),'omitnan');
                %%%%% selectivitydata
                rlmeans=[Selectivitydata.leftmean(k) Selectivitydata.rightmean(k)];
                Selectivitydata.score(k)=diff(rlmeans)/sum(rlmeans);

                % and selectivity index
                [~,Selectivitydata.pval(k)]=SelectivityIndex(spkevs(isCorr==2-k),odorid(isCorr==2-k),bootct); % 500 boots
                Selectivitydata.issig(k)=Selectivitydata.pval(k) <= odorPcrit & ~isnan(Selectivitydata.score(k)); % boolean
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
            catch
                fprintf('sess %d, cell %d didnt work for %s', ses, i, testing{k})
            end
        end
        SuperRat(ses).units(i).OdorRates=[spkevs, trialmat(:,3:4)];
        SuperRat(ses).units(i).OdorSelective=Selectivitydata;
        SuperRat(ses).units(i).startcurves=startcurves;
        SuperRat(ses).units(i).endcurves=endcurves;
    end
    fprintf('session %d done, took %.2f seconds \n', ses,toc)

end