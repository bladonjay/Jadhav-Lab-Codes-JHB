%% ClaireCellCellInteractions

%{
Theory: beta helps coordinate codes across PFC and CA1
Model: 
1. beta coherence increases suggesting a read frame for pyramidal cells
and that INs are responsible for synchronizing that read frame
 Prediction: to do this, INS and PYRs will have higher coherence to beta

2. This allows for better peformance: either
    a. better cue code
        Prediction: cue coding is stronger in correct trials durign high
        beta
    b. better translation from cue->response
    c. better response code
Previous datapoints:
1. Coherence in beta and rr band go up during odor sampling, but do not
change across correct and incorrect trials.  This suggests that while the
interneurons are oscillating coherently in each region, the information in
pyramidal cells may not be reached.
2. Cherence of INS in both CA1 and PFC are higher for correct trials,
coherence for PYR in CA1 is higher for correct trials


% some ideas:
1. see if spikes during odor and strong beta code better
2. see if spikes during odor and strong beta and right beta phase code
better
3. does the ensemble code shift during beta oscillations?

%}
%% previous work on timing
%{
Claire did a few plots that may be informative, lets see why these were
taken out of the analysis
1. CA_PFC beta coherence during odor presentation was anticorrelated with
sampling time: higher coherence==lower latency to choose

2. Neural discrimination preceded decision
    -CA1 coding occurred earlier

3. Neural population discrimination correlated with subsequent choice
    -I think the issue was that this was in the aggregate, and that one
    animal skewed the dataset because trials were averaged across


% so lets see if we can capture this beta timing analysis.
% some options here:
1. calculate on a trial to trial basis the first moment beta power is
significantly higher than chance. then calculate when the animal releases
his nose.  this should correlate on a trial to trial basis


% i dont think we need to worry too much about latency during this period,
I think correct vs incorrect may be the way to go here...  Need to see
whether cell xcorrs are different between correct and incorrect trials....
I think this could bear out as both cohere to beta better during correct
trials...
%}
%% okay so what claire does:

%{
1. gather all cells of your type
2. pull all spikes during the nosepoke window
3. run spiketrainxcorr
    a. bin spikes (count into 10 msec bins.
    b. from -250 to 250 msec each dir
    c. convolve w/gausswin of mean 20 msec (and size... 20 msec... this is wrong
    d. positive lags for her means cell 2 leads
4. after aggregating all cells, runs a one-sample ttest to see if the
peaks in each pairs lag are nonzero.
%}


% okay so the theory is....?
%{
1. beta oscillations help coordinate cross regional spiking.
    -beta coherence is high during odor time
    -beta locking is higher for correct trials
    -pyr-in interactions appear to be tight during odor period
    - does beta impact this?

options here:
- first replicate claires pyr-pyr and in-pyr analysis showing peaks in ac
- then figure out a way to measure whether this varies with beta, or
interacts with beta
- it mgiht naturally fall if cells are cohering to beta phases, btu how
would we show that?

%}

%% replicate the cell-cell interactions

% this turns the ts to a vector and xcorrs the vectors.  this wont take
% care of hanging zeroes... so dont use this method

%{

tbin=.005;
gauss=.01;
corrspan=.25;
span=-.25:tbin:.25;

minSample=.5;
celltypes={'pyr','in'};
regions={'PFC','CA1'};

plotIT=0;

allgrams=[1 1; 1 2; 2 1; 2 2]; varNames=[];
for i=1:length(allgrams); varNames{i}=[celltypes{allgrams(i,1)} '-' celltypes{allgrams(i,2)}]; end

corrTable=table('Size',[length(SuperRat),4],'VariableNames',varNames,'VariableTypes',repmat({'cell'},1,4));

for i=1:length(SuperRat)
    sessclock=tic;
    % save this out for each rat
    
    % only take correct trials
    
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    odorTimes=odorTimes(SuperRat(i).trialdata.CorrIncorr10==1,:);
    odorTimes(diff(odorTimes,1,2)>2.5 | diff(odorTimes,1,2)<minSample,:)=[]; % from 0.5 to 2 second samples
    odorTimes(:,1)=odorTimes(:,2)-minSample;
    % standardize trial length?
    % take last .5 seconds...
    
    
    % start with pyr-pyr and PFC is first region (as in regions
    for pr=1:length(varNames) % for each pairing
        PFCunits=SuperRat(i).units(cellfun(@(a) contains(a,regions{1}), {SuperRat(i).units.area}) & ...
            cellfun(@(a) contains(a,celltypes{allgrams(pr,1)}), {SuperRat(i).units.type}) & ...
            cellfun(@(a) any(a(:,4)), {SuperRat(i).units.OdorResponsive}));
        
        CA1units=SuperRat(i).units(cellfun(@(a) contains(a,regions{2}), {SuperRat(i).units.area}) & ...
            cellfun(@(a) contains(a,celltypes{allgrams(pr,2)}), {SuperRat(i).units.type}) & ...
            cellfun(@(a) any(a(:,4)), {SuperRat(i).units.OdorResponsive}));
        % corrmat will be big, it'll be a 3d, with i pfc cells by j HPC cells by k corrvals 
        corrMatrix=nan(length(PFCunits),length(CA1units),length(span)); % -25:25 bins
        
        if plotIT==1, figure; it=1; end 
        
        for j=1:length(PFCunits)
            Pspikes=event_spikes(PFCunits(j).ts(:,1), odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
            Pvec=SmoothMat2(accumarray(round(Pspikes/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);
            Pvecnull=SmoothMat2(accumarray(round(PFCunits(j).ts(:,1)/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);

            for k=1:length(CA1units)
                Cspikes=event_spikes(CA1units(k).ts, odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
                Cvec=SmoothMat2(accumarray(round(Cspikes/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);
                % nnull is all spikes
                Cvecnull=SmoothMat2(accumarray(round(CA1units(k).ts(:,1)/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);
                
                % i think for trialdata i may want to do this as square
                % matrices... lets see if we can do that
                corrMatrix(j,k,:)=xcorr(Pvec,Cvec,corrspan/tbin);
                
                if plotIT==1
                    subplot(length(PFCunits),length(CA1units),it);
                    plot(span,corrMatrix(j,k,:)); yyaxis right;
                    myxcorr=xcorr(Pvecnull,Cvecnull,corrspan/tbin);
                    plot(span,myxcorr);
                    it=it+1;
                end
            end
        end
        corrTable{i,pr}={corrMatrix};
        % on first blush, it looks like specific PFC units will assume
        % peaks in their xcorr with many CA1 units during odor sampling
        % she will plot the mean of all these
    end
    fprintf('Sess %d, rat %s %d done in %d, prob %d more\n',i,SuperRat(i).name,...
        SuperRat(i).daynum,round(toc(sessclock)), round(toc(sessclock)/i*(length(SuperRat)-i)));
end

%% now the workup of each celltype


for i=1:size(corrTable,2)
    allData=[];
    for j=1:size(corrTable,1)
        thisCorr=corrTable{j,i}{:};
        allData=[allData reshape(permute(thisCorr,[3,1,2]),size(thisCorr,3),size(thisCorr,1)*size(thisCorr,2))];
    end    
    allData=zscore(allData,1,1);
    figure('Position',[100+100*i 300 380 1000]);
    subplot(3,1,1);
    [a,b]=max(allData,[],1);
    [~,inds]=sort(b);
    imagesc(span,1:size(allData,2),allData(:,inds)');
     hold on; plot([0 0],[0 size(allData,2)],'k')
      title(sprintf('%s in PFC %s in CA1',celltype{allgrams(i,1)},celltype{allgrams(i,2)}));
    subplot(3,1,2);
    plot(span,mean(allData,2,'omitnan'));
    hold on;
    patch([span fliplr(span)], [mean(allData,2,'omitnan')+nanstd(allData,1,2); ...
        flipud([mean(allData,2,'omitnan')-nanstd(allData,1,2)])]','b',...
        'LineStyle','none','FaceAlpha',.3);
     hold on; plot([0 0],get(gca,'YLim'),'k')
    xlabel('PFC lag   PFC lead'); ylabel('Corr');
    xlim([span(1) span(end)]);
    subplot(3,1,3);
    [y,x]=histcounts(span(b),-corrspan:gauss/2:corrspan);
    plot(mean([x(1:end-1); x(2:end)]),SmoothMat2(y,[1 20],5));
     hold on; plot([0 0],get(gca,'YLim'),'k')
    xlim([-.15 .15]);
    title(sprintf('Signed-rank test p=%.2e',signrank(span(b))));
end
%}


%% replicate the cell-cell interactions with a square matrix cc to remake
% Claires Figure

% this uses the trial by trial xcorr and averages- not the best because
% trials with high spike counts will be weighted equally to those with
% few...


% Here I correct for edge effects.... by taking a trialwise correlation
tbin=.005;
gauss=.01;
corrspan=.25;
span=-corrspan:tbin:corrspan;
celltypes={'pyr','in'};
regions={'PFC','CA1'};

plotIT=0;

corrTable=table('Size',[length(SuperRat),4],'VariableNames',varNames,'VariableTypes',repmat({'cell'},1,4));
pairIndex=[1 1; 1 2; 2 1; 2 2]; varNames=[];
for i=1:length(pairIndex); varNames{i}=[celltypes{pairIndex(i,1)} '-' celltypes{pairIndex(i,2)}]; end
for i=1:length(SuperRat)
    sessclock=tic;
    % save this out for each rat
    

    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    odorTimes(diff(odorTimes,1,2)>2.5 | diff(odorTimes,1,2)<.5,:)=[]; % from 0.5 to 2 second samples
    % start with pyr-pyr and PFC is first region (as in regions
    for pr=1:length(varNames) % for each pairing
        PFCunits=find(cellfun(@(a) contains(a,regions{1}), {SuperRat(i).units.area}) & ...
            cellfun(@(a) contains(a,celltypes{pairIndex(pr,1)}), {SuperRat(i).units.type}) & ...
            cellfun(@(a) any(a(:,4)), {SuperRat(i).units.OdorResponsive}));
        
        CA1units=find(cellfun(@(a) contains(a,regions{2}), {SuperRat(i).units.area}) & ...
            cellfun(@(a) contains(a,celltypes{pairIndex(pr,2)}), {SuperRat(i).units.type}) & ...
            cellfun(@(a) any(a(:,4)), {SuperRat(i).units.OdorResponsive}));
        % corrmat will be big, it'll be a 3d, with i pfc cells by j HPC cells by k corrvals 
        odorMatrix=nan(length(PFCunits),length(CA1units),length(span)); % -25:25 bins
        
        if plotIT==1, figure; it=1; end 
        
        for j=1:length(PFCunits)
            [Pspikes,~,~,~,~,Pspk]=event_spikes(SuperRat(i).units(PFCunits(j)).ts(:,1), odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
            
            %Pvec=SmoothMat2(accumarray(round(Pspikes/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);

            %Pvecnull=SmoothMat2(accumarray(round(PFCunits(j).ts(:,1)/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);

            for k=1:length(CA1units)
                [Cspikes,~,~,~,~,Cspk]=event_spikes(SuperRat(i).units(CA1units(k)).ts, odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
                
                activetrials=cellfun(@(a,b) ~isempty(a) & ~isempty(b), Pspk, Cspk);
                Pvec=cellfun(@(a) SmoothMat2(accumarray(ceil(a/tbin),1),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    Pspk(activetrials), 'UniformOutput', false);
                Cvec=cellfun(@(a) SmoothMat2(accumarray(ceil(a/tbin),1),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    Cspk(activetrials), 'UniformOutput', false);
                %Cvec=SmoothMat2(accumarray(round(Cspikes/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);
                %Cvecnull=SmoothMat2(accumarray(round(CA1units(k).ts(:,1)/tbin),1),[gauss/tbin*10,1],[gauss/tbin,1]);
                
                % i think for trialdata i may want to do this as square
                % matrices... lets see if we can do that
                % trials ^2
                it=1;
                trCorr=[];
                for xt=1:length(Pvec)
                    trCorr(it,:)=xcorr(Pvec{xt},Cvec{xt},corrspan/tbin);
                    it=it+1;
                end
                % whats the null?
                % null is boot circshift same trials
                
                %corrMatrix(j,k,:)=xcorr(Pvec,Cvec,corrspan/tbin);
                odorMatrix(j,k,:)=mean(zscore(trCorr,[],2),1,'omitnan');
                
                if plotIT==1
                    subplot(length(PFCunits),length(CA1units),it);
                    plot(span,odorMatrix(j,k,:)); yyaxis right;
                    myxcorr=xcorr(Pvecnull,Cvecnull,corrspan/tbin);
                    plot(span,myxcorr);
                    it=it+1;
                end
            end
        end
        corrTable{i,pr}={odorMatrix};
        % on first blush, it looks like specific PFC units will assume
        % peaks in their xcorr with many CA1 units during odor sampling
        % she will plot the mean of all these
    end
    fprintf('Sess %d, rat %s %d done in %d, prob %d more\n',i,SuperRat(i).name,...
        SuperRat(i).daynum,round(toc(sessclock)), round(toc(sessclock)/i*(length(SuperRat)-i)));
end

%% now the workup of each celltype


for i=1:size(corrTable,2)
    allData=[];
    for j=1:size(corrTable,1)
        thisCorr=corrTable{j,i}{:};
        allData=[allData reshape(permute(thisCorr,[3,1,2]),size(thisCorr,3),size(thisCorr,1)*size(thisCorr,2))];
    end    
    allData=zscore(allData,1,1);
    figure('Position',[100+100*i 300 380 1000]);
    subplot(3,1,1);
    [a,b]=max(allData,[],1);
    [~,inds]=sort(b);
    imagesc(span,1:size(allData,2),allData(:,inds)');
     hold on; plot([0 0],[0 size(allData,2)],'k')
      title(sprintf('%s in PFC %s in CA1',celltype{pairIndex(i,1)},celltype{pairIndex(i,2)}));
    subplot(3,1,2);
    plot(span,mean(allData,2,'omitnan'));
    hold on;
    patch([span fliplr(span)], [mean(allData,2,'omitnan')+nanstd(allData,1,2); ...
        flipud([mean(allData,2,'omitnan')-nanstd(allData,1,2)])]','b',...
        'LineStyle','none','FaceAlpha',.3);
     hold on; plot([0 0],get(gca,'YLim'),'k')
    xlabel('PFC lag   PFC lead'); ylabel('Corr');
    xlim([span(1) span(end)]);
    subplot(3,1,3);
    [y,x]=histcounts(span(b),-corrspan:gauss/2:corrspan);
    plot(mean([x(1:end-1); x(2:end)]),SmoothMat2(y,[1 20],5));
     hold on; plot([0 0],get(gca,'YLim'),'k')
    xlim([-.15 .15]);
    title(sprintf('Signed-rank test p=%.2e',signrank(span(b))));
end

%%
%
%
%
%
%    Identify cell pairs who have significant interactions
%
%
%
%
%
%{
The theoretical question is whether cell cell interactions between PFC and
CA1 show evidence for cross regional information flow during odor sampling
and whether this is different for correct and incorrect trials

the initial suggestion is that the cells fail to interact and therefore
information is not passed between regions effectively

One alternative is that ca1->PFC interactions fail to take place

a second alternative is that PFC->CA1 interactions fail to take place

a third alternative is that during fails, there is too much top=down and
not enough bottom up

 to analyze these- we need to
1. determine significantly connected pairs
2. get the direction of that connectivity

3. for a connected pair, does that direction change? say from sampling to
running?

%}


%% xcorr for odor responsive cell pairs with null

% Claires Figure
% this is for the odor period
% this also works up all odor responsive cells
% the question being what do the xcorrgrams look like for responsive cells
% vs selective cells... can we see evidence of beta in selective cells, or
% do we see peaks between selective cells only?

% Here I correct for edge effects.... by taking a trialwise correlation
tbin=.005;
gauss=.015;
corrspan=.25;
span=-corrspan:tbin:corrspan;
regions={'PFC','CA1'};
nboots=500;
plotIT=0;
odorDur=0.5;
runDur=4;

corrTable=table('Size',[length(SuperRat),4],'VariableNames',varNames,'VariableTypes',repmat({'cell'},1,4));
pairIndex=[1 1; 1 2; 2 1; 2 2]; varNames=[];
for i=1:length(pairIndex); varNames{i}=[celltypes{pairIndex(i,1)} '-' celltypes{pairIndex(i,2)}]; end


for i=1:length(SuperRat)
    sessclock=tic;
    % save this out for each rat
    
    
    % start with pyr-pyr and PFC is first region (as in regions
    PFCunits=find(cellfun(@(a) contains(a,regions{1}), {SuperRat(i).units.area}) & ...
        cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
    
    CA1units=find(cellfun(@(a) contains(a,regions{2}), {SuperRat(i).units.area}) & ...
        cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
    
    for q=[PFCunits CA1units]
        SuperRat(i).units(q).partners=[];
    end
    
    
    % all trials, not c vs i yet
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    activetrials=diff(odorTimes,1,2)<2.5 & diff(odorTimes,1,2)>odorDur; % from 0.5 to 2 second samples
    
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10; % from 0.5 to 2 second samples
    
    % and get runtimes
    try
        runTimes=[SuperRat(i).trialdata.sniffend+odorDur SuperRat(i).trialdata.rewardstart];
    catch
        try
            runTimes=[SuperRat(i).trialdata.sniffend+odorDur SuperRat(i).trialdata.endtime];
        catch
            continue
        end
    end
    
    % use only okay trials
    odorTimes=odorTimes(activetrials,:);
    trialCorrect=trialCorrect(activetrials,:);
    runTimes=runTimes(activetrials,:);
    
    
    % cast all odorTimes to 0.5 seconds
    odorTimes(:,1)=odorTimes(:,2)-odorDur;
    % for run go forward 4 seconds
    runTimes(:,2)=runTimes(:,1)+runDur;
    
    
    
    
    % corrmat will be big, it'll be a 3d, with i pfc cells by j HPC cells by k corrvals
    odorMatrix=nan(length(PFCunits),length(CA1units),length(span)); % -25:25 bins
    
    for j=1:length(PFCunits)
        [Pspikes,~,~,~,~,Pspk]=event_spikes(SuperRat(i).units(PFCunits(j)).ts(:,1), odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
        %[Prspikes,~,~,~,~,Prspk]=event_spikes(SuperRat(i).units(PFCunits(j)).ts(:,1), runTimes(:,1),0,runTimes(:,2)-runTimes(:,1));
        
        
        for k=1:length(CA1units)
            [Cspikes,~,~,~,~,Cspk]=event_spikes(SuperRat(i).units(CA1units(k)).ts, odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
            %[Crspikes,~,~,~,~,Crspk]=event_spikes(SuperRat(i).units(CA1units(k)).ts(:,1), runTimes(:,1),0,runTimes(:,2)-runTimes(:,1));
            
            oktrials=cellfun(@(a,b) ~isempty(a) & ~isempty(b), Pspk, Cspk);
            
            if sum(oktrials)>8
                % need to do odor and run, and randomize each, and then do
                % correct vs incorrect and jackknife an effect size
                % need to make the vectors the same length
                Pvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(.5/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    Pspk(oktrials), 'UniformOutput', false);
                %Prvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(4/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                %    Prspk(activetrials2), 'UniformOutput', false);
                
                Cvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(.5/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    Cspk(oktrials), 'UniformOutput', false);
                %Crvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(4/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                %    Crspk(activetrials2), 'UniformOutput', false);
                
                % i think for trialdata i may want to do this as square
                % matrices... lets see if we can do that
                % neg values means second vector is shfited in the -
                % direction, and so high corr means top leads bottom
                % neg values means PFC leads
                it=1;
                trCorr=[]; trCorrR=[];
                
                % minus values means that the first vector (PFC) leads
                for xt=1:length(Pvec)
                    trCorr(it,:)=xcorr(Pvec{xt},Cvec{xt},corrspan/tbin);
                    %trCorrR(it,:)=xcorr(Prvec{xt},Crvec{xt},corrspan/tbin);
                    it=it+1;
                end
                % for boot
                bootCorr=[]; %bootCorrR=[];
                for bt=1:nboots
                    mixup=Pvec(randperm(length(Pvec)));
                    %mixup2=Prvec(randperm(length(Prvec)));
                    %randshift=randi(round(corrspan/tbin*2*[.1 .9]),length(Pvec),1);
                    trCorrBtR=[]; trCorrBt=[];
                    for xt=1:length(Pvec)
                        trCorrBt(xt,:)=xcorr(mixup{xt},Cvec{xt},corrspan/tbin);
                        %trCorrBt(xt,:)=xcorr(circshift(Pvec{xt},randshift(xt))',Cvec{xt},corrspan/tbin);
                        %trCorrBtR(xt,:)=xcorr(mixup2{xt},Crvec{xt},corrspan/tbin);
                        
                    end
                    bootCorr(bt,:)=mean(trCorrBt,'omitnan');
                    %bootCorrR(bt,:)=mean(trCorrBtR,'omitnan');
                end
                
                % now find best peak within +/- 50 msec and a z score above
                % 2
                bootZ=(mean(trCorr,'omitnan')-mean(bootCorr,'omitnan'))./std(bootCorr,[],1,'omitnan');
                [mypeaks]=findpeaks(bootZ,2); mypeaks.loc=mypeaks.loc(abs(span(mypeaks.loc))<.1);
                
                % if theres an index where z>2
                if ~isempty(mypeaks.loc)
                    peakoffset=span(min(mypeaks.loc));
                else
                    peakoffset=nan;
                end
                PFCid=[SuperRat(i).units(PFCunits(j)).tet SuperRat(i).units(PFCunits(j)).unitnum];
                PFCtype=SuperRat(i).units(PFCunits(j)).type;
                PFCsel=SuperRat(i).units(PFCunits(j)).OdorSelective; % odor selective here
                PFCbeta=SuperRat(i).units(PFCunits(j)).betaRstat(1);
                PFCrr=SuperRat(i).units(PFCunits(j)).respRstat(1);
                
                CA1id=[SuperRat(i).units(PFCunits(j)).tet SuperRat(i).units(CA1units(k)).unitnum];
                CA1type=SuperRat(i).units(CA1units(k)).type;
                CA1sel=SuperRat(i).units(CA1units(k)).OdorSelective;
                CA1beta=SuperRat(i).units(CA1units(k)).betaRstat(2);
                CA1rr=SuperRat(i).units(CA1units(k)).respRstat(2);
               
                % now add to each cells partners
                SuperRat(i).units(PFCunits(j)).partners=[SuperRat(i).units(PFCunits(j)).partners; ...
                    struct('offset',peakoffset,'ID',CA1id,'selective',CA1sel,'type',CA1type,'gram',bootZ,...
                    'Bcoh',CA1beta,'RRcoh',CA1rr)];
                SuperRat(i).units(CA1units(k)).partners=[SuperRat(i).units(CA1units(k)).partners; ...
                    struct('offset',peakoffset,'ID',PFCid,'selective',PFCsel,'type',PFCtype,...
                    'Bcoh',PFCbeta,'RRcoh',PFCrr)];
                if ~isnan(peakoffset)
                    fprintf('Found one PFC%s to CA1%s pair, offset %.f msec \n',SuperRat(i).units(PFCunits(j)).type,...
                        SuperRat(i).units(CA1units(k)).type,peakoffset*1000);
                end
                
                if plotIT==1
                    figure;
                    subplot(2,2,1);
                    plot(span,mean(trCorr,'omitnan'));
                    hold on;
                    patch([span,fliplr(span)], [prctile(bootCorr,99), ...
                        fliplr(prctile(bootCorr,1))],[.4 .4 .4],...
                        'LineStyle','none','FaceAlpha',.4);
                    title(sprintf('CA1 %s %d PFC %s %d',SuperRat(i).units(CA1units(k)).type, CA1units(k),...
                        SuperRat(i).units(PFCunits(j)).type, PFCunits(j)));
                    xlabel('PFC leads CA1 leads')
                    subplot(2,2,2);
                    plot(span,mean(trCorrR,'omitnan'));
                    hold on;
                    patch([span,fliplr(span)], [mean(bootCorrR,'omitnan')+nanstd(bootCorrR), ...
                        fliplr(mean(bootCorrR,'omitnan')-nanstd(bootCorrR))],[.4 .4 .4],...
                        'LineStyle','none','FaceAlpha',.4);
                    title(sprintf('CA1 %s %d PFC %s %d',SuperRat(i).units(CA1units(k)).type, CA1units(k),...
                        SuperRat(i).units(PFCunits(j)).type, PFCunits(j)));
                    xlabel('PFC leads CA1 leads')
                    subplot(2,2,3);
                    plot(span,(mean(trCorr,'omitnan')-mean(bootCorr,'omitnan'))./std(bootCorr,[],1,'omitnan'));
                    title('odor period');
                    subplot(2,2,4);
                    plot(span,(mean(trCorrR,'omitnan')-mean(bootCorrR,'omitnan'))./std(bootCorrR,[],1,'omitnan'));
                    title('Run period');
                end
            end
        end
    end
    corrTable{i,pr}={odorMatrix};
    % on first blush, it looks like specific PFC units will assume
    % peaks in their xcorr with many CA1 units during odor sampling
    % she will plot the mean of all these
    fprintf('Sess %d, rat %s %d done in %d, prob %d more\n',i,SuperRat(i).name,...
        SuperRat(i).daynum,round(toc(sessclock)), round(toc(sessclock)/i*(length(SuperRat)-i)));
end



%% 

% simple question is this... is there any relationship between coding, cell
% pairs and oscillations.

% so to gather:
%{
1. coding properties
2. cell partners
3. oscillation coherence
%}


% Okay, I see no interaction between these features... lets make sure the
% xcorr is valid though...
%% now we count the overlaps
    
    % table will be: tot, beta coh, rr coh, selective, hasPartners
    varnames={'total','beta','resp','bothR','selective','selectiveB','selectiveRR',...
        'hasPartners','PartnerSelective','PartnerB','PartnerR'};
% first n partners vs coherence
    barTable=table('Size',[2,11],'VariableTypes',repmat({'double'},11,1),...
        'VariableNames',varnames);
    
hf1=figure; 
hf2=figure;
regions={'PFC','CA1'};
celltype={'pyr','in'}; % will need to breakout celltypes soon
ct=1;
for reg=1:length(regions)
    % PFC first
    megarat=[];
    
    for i=1:length(SuperRat)
        theseUnits=find(cellfun(@(a) contains(a,regions{reg}), {SuperRat(i).units.area}) & ...
            cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive})); % &...
            %cellfun(@(a) contains(a,celltype{ct}), {SuperRat(i).units.type}));
        if ~isempty(theseUnits)
        % gather n partners,
        if isfield(SuperRat(i).units,'RUNpartners')
        megarat=[megarat rmfield(SuperRat(i).units(theseUnits),'RUNpartners')];
        else
            megarat=[megarat SuperRat(i).units(theseUnits)];
        end
        end
    end
    %megarat(cellfun(@(a) isempty(a), {megarat.partners}))=[];
    % now tabulate the parameters
    nPartners=[];
    for i=1:length(megarat)
        try
         nPartners(i)=sum(~isnan(cell2mat({megarat(i).partners.offset})));
        catch
            nPartners(i)=0;
        end
    end
    hasPeer=nPartners>0;
    %offsets=cellfun(@(a) a, {megarat.partners});
    betaMVL=cellfun(@(a) a(reg), {megarat.betaMVL});
    betaP=cellfun(@(a) a(reg)<.05, {megarat.betaRstat});
    respMVL=cellfun(@(a) a(reg), {megarat.respMVL});
    respP=cellfun(@(a) a(reg)<.05, {megarat.respRstat});
    mySel=cellfun(@(a) a(2), {megarat.OdorSelective});
    isSel=cellfun(@(a) a(4)<.05, {megarat.OdorSelective});
    

    figure(hf1);
    subplot(2,3,(reg-1)*3+1);
    % first we do beta and rr

    barTable.total(reg)=length(nPartners);
    barTable.beta(reg)=sum(betaP);
    barTable.resp(reg)=sum(respP);
    barTable.bothR(reg)=sum(betaP & respP);
    
    props=barTable{reg,2:4}/barTable{reg,1};
    b=bar(1:4,[props props(1)*props(2)],'LineStyle','none','FaceColor','flat');
    [~,bounds]=binofit(ceil(props(2)*props(3)*barTable.total(reg)),barTable.total(reg));
    b.CData(4,:)=[.7 .7 .7];
    hold on; errorbar(4, props(1)*props(2), bounds(1)'-props(1)*props(2),bounds(2)'-props(1)*props(2),'k.','CapSize',0,'LineWidth',2)
    set(gca,'XTick',[1:4],'XTickLabel',{'Beta','RR','both','both null'});
    title(sprintf(' %s Cells dont cohere to both beta and RR',regions{reg}));
    box off
    
    barTable.selective(reg)=sum(isSel);
    barTable.selectiveB(reg)=sum(betaP & isSel);
    barTable.selectiveRR(reg)=sum(respP & isSel);
    
    
    subplot(2,3,(reg-1)*3+2);
    props=barTable{reg,5:7}/barTable{reg,1};
    nulls=barTable{reg,[2 3]}/barTable{reg,1}*(barTable{reg,[5]}/barTable{reg,1});
    b=bar(1:5,[props(1:2) nulls(1) props(3) nulls(2)],'LineStyle','none','FaceColor','flat');
    [~,bounds]=binofit(ceil(nulls*barTable.total(reg)),barTable.total(reg));
    b.CData([3 5],:)=[.7 .7 .7; .7 .7 .7];

    hold on; errorbar([3 5], nulls, bounds(:,1)'-nulls,bounds(:,2)'-nulls,'k.','CapSize',0,'LineWidth',2)
    set(gca,'XTick',[1:5],'XTickLabel',{'Selective','Beta & Sel','B & sel chance','RR & Sel','RR & sel chance'});
    box off;
    title('Selective Cells dont Cohere More');
    
    barTable.hasPartners(reg)=sum(nPartners>0);
    barTable.PartnerSelective(reg)=sum(nPartners>0 & isSel);
    barTable.PartnerB(reg)=sum(nPartners>0 & betaP);
    barTable.PartnerR(reg)=sum(nPartners>0 & respP);
    
    
    subplot(2,3,(reg-1)*3+3);
    props=barTable{reg,8:11}/barTable{reg,1};
    nulls=props(1)*(barTable{reg,[5 2 3]}/barTable{reg,1});
    b=bar([1 1.1 2:7],linearize([props; nan nulls])',6,'LineStyle','none','FaceColor','flat');
    b.CData([4 6 8],:)=repmat([.7 .7 .7],3,1);

    [~,bounds]=binofit(ceil(nulls*barTable.total(reg)),barTable.total(reg));
    hold on; errorbar([3 5 7], nulls, bounds(:,1)'-nulls,bounds(:,2)'-nulls,'k.','CapSize',0,'LineWidth',2)
    set(gca,'XTick',[1:7],'XTickLabel',{'Peers','Peers & Sel','Peers & Sel Null','Peers & Beta',...
        'Peers & Beta Null','Peers & RR','Peers & RR Null'});
    box off;
    title('cells w partners might be more selective & rr coherent?');

    %{
    % if you want to look at conditional probabilities...
    figure(hf2);
 
    % now lets look at conditional probabilities
    % first resp cells, %% of those that are beta too
    
    %mybars=[sum(~respP(betaP)) sum(~betaP(respP)); sum(respP(betaP)) sum(betaP(respP))]/length(betaP);
    subplot(1,2,reg); bar(mybars','stacked')
    set(gca,'XTickLabel',{'beta Coherent','RR Coherent'});
    legend('','both rhythms');
    title(regions{reg});
    
    % prop of coherent cells are selective
    mybars=[sum(~isSel(betaP)) sum(~isSel(respP)) sum(isSel); sum(isSel(betaP)) sum(isSel(respP)) 0]/length(betaP);
    subplot(1,2,reg); bar(mybars','stacked')
    set(gca,'XTickLabel',{'beta Coherent','RR Coherent'});
    legend('nonsel','selective');
    title(regions{reg});
    
    % prop of selective cells are coherent
    mybars=[sum(~betaP(isSel) & ~respP(isSel)) sum(betaP) sum(respP);...
       sum(betaP(isSel) & ~respP(isSel)) 0 0;...
        sum(~betaP(isSel) & respP(isSel)) 0 0;...
        sum(betaP(isSel) & respP(isSel)) 0 0]/length(betaP);
    subplot(1,2,reg); bar(mybars','stacked')
    set(gca,'XTickLabel',{'Selective','Bcoh','RRcoh'});
    legend('not coherent','only beta', 'only rr','both');
    title(regions{reg});
    linkaxes(get(gcf,'Children'),'y')
    
    % prop of cells w peers are coherent
    mybars=[sum(~betaP(hasPeer) & ~respP(hasPeer)) sum(betaP) sum(respP);...
        sum(betaP(hasPeer) & ~respP(hasPeer)) 0 0;...
         sum(~betaP(hasPeer) & respP(hasPeer)) 0 0;...
         sum(betaP(hasPeer) & respP(hasPeer)) 0 0]/length(betaP);
    subplot(1,2,reg); bar(mybars','stacked')
    set(gca,'XTickLabel',{'hasPeers','Bcoh','RRcoh'});
    legend('not coherent','only beta', 'only rr','both');
    title(regions{reg});
    linkaxes(get(gcf,'Children'),'y')
    

   % prop of peer cells are selective
    mybars=[sum(~hasPeer(isSel)) sum(~isSel(hasPeer)); sum(hasPeer(isSel)) sum(isSel(hasPeer)) ]/length(betaP);
    subplot(1,2,reg); bar(mybars','stacked')
    set(gca,'XTickLabel',{'Selective','hasPeers'});
    legend('Only', 'has both');
    title(regions{reg});
    ylabel('Proportion of Task Responsive Units')
%}
end
    
linkaxes(get(hf2,'Children'),'y')

%% do cells pairs that independently cohere to beta or rr have different offsets?



% get PFC cells and go PFC->HPC becaue thats the direction of the ccgram


    
    
%%
%
% Lets get some xcorrellogram characteristics
%
%
%
%
%
%
%
%
%%
%do a fourier on the spike xcorr, get a beta to theta ratio

% then plot that against the direction of the center.

%ken kay removed the slow componenet to get the center and the xcorr

% get the offset by taking the highest point within +/-50 msec

%% Gather the cell-cell xcorrellograms, find crucial stats

% matt jones and matt wilson, 1005 plos theta rhythms coordinate hc pfc
% interactions in spatial memory task

% basically they measure the peak from +200msec to -200 msec and do a
% difference score

% I also want to measure the amount of beta or theta in the xcgram for c
% and incorrect

% use the get ACinfo to measure beta and theta in autocorr, do this for
% correct first just to see what most pairs do, then see if it's different
% for incorrect trials (will need to do some stats here)


% Here I correct for edge effects.... by taking a trialwise correlation
tbin=.005;  gauss=.015;  corrspan=.4;
span=-corrspan:tbin:corrspan;
regions={'PFC','CA1'};
nboots=200;  plotIT=0;  odorDur=0.5;

pairIndex=[0 0; 1 0; 0 1; 1 1];

varNames={'xcorrgrams','statstable','celltypes'};

corrTable=table('Size',[length(SuperRat),3],'VariableNames',varNames,'VariableTypes',repmat({'cell'},1,3));

for i=1:length(SuperRat)
    sessclock=tic;
       
    % start with all cells and PFC is first region (as in regions)
    PFCunits=find(cellfun(@(a) contains(a,regions{1}), {SuperRat(i).units.area}) & ...
        cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
    
    CA1units=find(cellfun(@(a) contains(a,regions{2}), {SuperRat(i).units.area}) & ...
        cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
    
    
    
    % all trials, not c vs i yet
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    oktrials=diff(odorTimes,1,2)<2.5 & diff(odorTimes,1,2)>.5; % from 0.5 to 2 second samples
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10; % from 0.5 to 2 second samples
    
    % use only okay trials
    odorTimes=odorTimes(oktrials,:);
    trialCorrect=trialCorrect(oktrials,:);
    
    % cast all odorTimes to 0.5 seconds
    odorTimes(:,1)=odorTimes(:,2)-odorDur;
    % corrmat will be big, it'll be a 3d, with i pfc cells by j HPC cells by k corrvals
    odorMatrix=nan(length(PFCunits),length(CA1units),length(span)); % -25:25 bins
    corrMatrixX=odorMatrix;
    odorStats=nan(length(PFCunits),length(CA1units),3); % this will be center in msec, and theta-gamma ratio

    cellTypes=nan(length(PFCunits),length(CA1units));
    % celltypetable- 0=IN IN; 1=PYR-IN; 2=IN-PYR; 3=PYR-PYR
    for j=1:length(PFCunits)
        [Pspikes,~,~,~,~,Pspk]=event_spikes(SuperRat(i).units(PFCunits(j)).ts(:,1),...
            odorTimes(trialCorrect==1,1),0,odorTimes(trialCorrect==1,2)-odorTimes(trialCorrect==1,1));
        [PXspikes,~,~,~,~,PXspk]=event_spikes(SuperRat(i).units(PFCunits(j)).ts(:,1),...
            odorTimes(trialCorrect==0,1),0,odorTimes(trialCorrect==0,2)-odorTimes(trialCorrect==0,1));
        Ptype=contains(SuperRat(i).units(PFCunits(j)).type,'pyr');
        
        for k=1:length(CA1units)
            [Cspikes,~,~,~,~,Cspk]=event_spikes(SuperRat(i).units(CA1units(k)).ts,...
                odorTimes(trialCorrect==1,1),0,odorTimes(trialCorrect==1,2)-odorTimes(trialCorrect==1,1));
            [CXspikes,~,~,~,~,CXspk]=event_spikes(SuperRat(i).units(CA1units(k)).ts,...
                odorTimes(trialCorrect==0,1),0,odorTimes(trialCorrect==0,2)-odorTimes(trialCorrect==0,1));
            
            Ctype=contains(SuperRat(i).units(CA1units(k)).type,'pyr');
            
            cellTypes(j,k)=find(sum([Ptype Ctype]==pairIndex,2)==2)-1;

            
            activetrials=cellfun(@(a,b) ~isempty(a) & ~isempty(b), Pspk, Cspk);
            activetrials2=cellfun(@(a,b) ~isempty(a) & ~isempty(b), PXspk, CXspk);
            
            if sum(activetrials)>=4
                % need to do Correct and Incorrect, and randomize each, and then do

                % i smooth before, claire smoothes after
                Pvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(odorDur/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    Pspk(activetrials), 'UniformOutput', false);
                PXvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(odorDur/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    PXspk(activetrials2), 'UniformOutput', false);
                
                Cvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(odorDur/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    Cspk(activetrials), 'UniformOutput', false);
                CXvec=cellfun(@(a) SmoothMat2(histcounts(ceil(a/tbin),[0:ceil(odorDur/tbin)]'),[1,gauss/tbin*10],[1,gauss/tbin]),...
                    CXspk(activetrials2), 'UniformOutput', false);
                
                % neg values means PFC leads
                it=1;
                trCorr=[]; trCorrX=[];
                for xt=1:length(Pvec)
                    trCorr(it,:)=xcorr(Pvec{xt},Cvec{xt},corrspan/tbin);
                    it=it+1;
                end
                for xt=1:length(PXvec)
                    trCorrX(it,:)=xcorr(PXvec{xt},CXvec{xt},corrspan/tbin);
                    it=it+1;
                end
                
                
                % for boot
                bootCorr=[]; bootCorrX=[];

                for bt=1:nboots
                    mixup=Pvec(randperm(length(Pvec)));
                    mixup2=PXvec(randperm(length(PXvec)));
                    trCorrBtX=[]; trCorrBt=[];
                    for xt=1:length(Pvec)
                        trCorrBt(xt,:)=xcorr(mixup{xt},Cvec{xt},corrspan/tbin);
                    end
                    for xt=1:length(PXvec)
                        trCorrBtX(xt,:)=xcorr(mixup2{xt},CXvec{xt},corrspan/tbin);
                    end
                    bootCorr(bt,:)=mean(trCorrBt,'omitnan');
                    bootCorrX(bt,:)=mean(trCorrBtX,'omitnan');
                end

                %normcurve=xcorr(cell2mat(Pvec),cell2mat(Cvec),corrspan/tbin);
                normcurve=rescale(mean(trCorr,'omitnan')');
                %{
                mymodel= @(a,b,c,w,tau1,tau2,off,x)...
                    (a*(cos(2*pi*w*(x+off))+1)+b).*exp(-abs((x+off))/tau1) + c*exp(-((x+off).^2)/tau2^2);
                [mx] = fit(span', normcurve, mymodel, ...
                    'StartPoint', [.5,  .1,   3,   7,  .2, .01, 0], ...
                    'Lower',      [ 0,  0,   0,   4,  .1,  0, -.08], ...
                    'Upper',      [10,  10,  10, 12,  10, .1, .08],...
                    'MaxFunEvals',10^10,'TolFun',10^-10, 'Robust', 'LAR');
                
                [mx2] = fit(span', normcurve, mymodel, ...
                    'StartPoint', [ .5,  .1,   3,  22,  .2, .01,   0], ...
                    'Lower',      [  0,   0,   0,  20,  .1,   0, -.08], ...
                    'Upper',      [ 10,  10,  10,  35,   10,  .1, .08],...
                    'MaxFunEvals',10^10,'TolFun',10^-10, 'Robust', 'LAR');
                
                cf=coeffvalues(mx); cf2=coeffvalues(mx2);
                
                corrStats(j,k,:)=[nan cf(7) (cf(1)-cf(2))/(cf2(1)+cf2(2))];
                %}
                % or run fourier analysis on the acgram
                %params=struct('tapers',[3 5],'Fs',1/tbin,'pad',1,'fpass',[2 40]);
                % need at least 0.4 sec of sample to get 5 hz
                [pxx,f] = pwelch(normcurve,round(.4/tbin),round(.3/tbin),5:.5:50,1/tbin,'Power');
                %[S,f]=mtspectrumc(normcurve',params)
                %[pxx,f] = pmtm(normcurve,4,length(normcurve),1/tbin);
                spect=rescale(pxx.*f);
                thetapwr=max(spect(f>6 & f<12));
                betapwr=max(spect(f>20 & f<35));

                % now find best peak within +/- 100 msec
                bootZ=(mean(trCorr,'omitnan')-mean(bootCorr,'omitnan'))./std(bootCorr,[],1,'omitnan');
                [mypeaks]=findpeaks(abs(bootZ),2); mypeaks.loc(abs(span(mypeaks.loc))>.1)=[];
                [~,ind]=min(abs(span(mypeaks.loc))); % pick best peak
                
                if isempty(ind) 
                     odorStats(j,k,:)=[nan nan nan];
                else
                    % z value, time offset, and the beta/theta power ratio-
                    % theta pwr large, 3rd val is above zero, beta pwr
                    % large, 3rd val below zero
                    odorStats(j,k,:)=[bootZ(mypeaks.loc(ind)) span(mypeaks.loc(ind)) (thetapwr-betapwr)/(thetapwr+betapwr)];
                    %corrStats(j,k,:)=[cf(3) cf(7) (cf(1)-cf(2))/(cf2(1)+cf2(2))];
                end
                
                odorMatrix(j,k,:)=mean(trCorr,1,'omitnan');
                %corrMatrix(j,k,:)=bootZ;

                corrMatrixX(j,k,:)=mean(trCorr,1,'omitnan');
                
                if ~isempty(ind) % if there is a peak, log it
                    
                    
                    if plotIT==1
                        %{
                        figure;
                        subplot(2,2,1);
                        plot(span,mean(trCorr,'omitnan'));
                        hold on;
                        patch([span,fliplr(span)], [prctile(bootCorr,99), ...
                            fliplr(prctile(bootCorr,1))],[.4 .4 .4],...
                            'LineStyle','none','FaceAlpha',.4);
                        title(sprintf('CA1 %s %d PFC %s %d',SuperRat(i).units(CA1units(k)).type, CA1units(k),...
                            SuperRat(i).units(PFCunits(j)).type, PFCunits(j)));
                        xlabel('PFC leads CA1 leads')
                        subplot(2,2,2);
                        plot(span,mean(trCorrX,'omitnan'));
                        hold on;
                        patch([span,fliplr(span)], [prctile(bootCorrX,99), ...
                            fliplr(prctile(bootCorrX,1))],[.4 .4 .4],...
                            'LineStyle','none','FaceAlpha',.4);
                        title(sprintf('CA1 %s %d PFC %s %d',SuperRat(i).units(CA1units(k)).type, CA1units(k),...
                            SuperRat(i).units(PFCunits(j)).type, PFCunits(j)));
                        xlabel('PFC leads CA1 leads')
                        subplot(2,2,3);
                        plot(span,(mean(trCorr,'omitnan')-mean(bootCorr,'omitnan'))./std(bootCorr,[],1,'omitnan'));
                        title('odor period C');
                        subplot(2,2,4);
                        plot(span,(mean(trCorrX,'omitnan')-mean(bootCorrX,'omitnan'))./std(bootCorrX,[],1,'omitnan'));
                        title('odor period I');
                        %}
                        
                        
                        figure; subplot(2,1,1); 
                        plot(span,normcurve); yyaxis right;
                        plot(span,bootZ);
                        %{
                        subplot(4,1,2);
                        slow=(cf(1)*(cos(2*pi*cf(4)*(x+cf(7)))+1)+cf(2)).*exp(-abs((x+cf(7)))/cf(5));
                        fast=cf(3)*exp(-((x+cf(7)).^2)/cf(6)^2);
                        
                        b=bar(span,normcurve,'FaceColor','flat','EdgeColor','flat');
                        b.CData=[.7 .7 .7]; hold on;
                        plot(x,slow,'k'); plot(x,fast,'r'); plot(x,slow+fast,'b');
                        legend('data','slowRR','fast','fit');
                          slow=(cf2(1)*(cos(2*pi*cf2(4)*(x+cf2(7)))+1)+cf2(2)).*exp(-abs((x+cf2(7)))/cf2(5));
                        fast=cf2(3)*exp(-((x+cf2(7)).^2)/cf2(6)^2);
                        subplot(4,1,3);

                        b=bar(span,normcurve,'FaceColor','flat','EdgeColor','flat');
                        b.CData=[.7 .7 .7]; hold on;
                        plot(x,slow,'k'); plot(x,fast,'r'); plot(x,slow+fast,'b');
                        legend('data','slowRR','fast','fit');
                        %}
                        subplot(2,1,2);plot(f(f>5),spect(f>5));
                        title(sprintf('RR = %.2e, beta = %.2e',thetapwr,betapwr));
                    end
                end
            end
        end
    end
    
    corrTable.xcorrgrams(i)={odorMatrix};
    corrTable.corrStats(i)={odorStats};
    corrTable.cellTypes(i)={cellTypes};
    % on first blush, it looks like specific PFC units will assume
    % peaks in their xcorr with many CA1 units during odor sampling
    % she will plot the mean of all these
    fprintf('Sess %d, rat %s %d done in %d, prob %d more\n',i,SuperRat(i).name,...
    SuperRat(i).daynum,round(toc(sessclock)), round(toc(sessclock)/i*(length(SuperRat)-i)));
end
%%
% first want to produce waterfall plots separately for pairs that have
% theta, and cells that have beta.


% then I want to do some statistics on the offsets of their peaks... also
% wonder which have significant peaks

% first is PFC, second is CA1
% celltypetable- 0=IN IN; 1=PYR-IN; 2=IN-PYR; 3=PYR-PYR

corrTable(cellfun(@(a) isempty(a), corrTable.xcorrgrams),:)=[];
pairIndex=[]; allstats=[]; allpairs=[];
for i=1:height(corrTable)
    allGrams=corrTable.xcorrgrams{i};
    thisStats=corrTable.corrStats{i};
    pairTypes=corrTable.cellTypes{i};
    allstats=[allstats reshape(permute(thisStats,[3,1,2]),...
        [3,size(thisStats,1)*size(thisStats,2)])];
    pairIndex=[pairIndex reshape(permute(allGrams,[3,1,2]),...
        [size(corrTable.xcorrgrams{1},3),size(allGrams,1)*size(allGrams,2)])];
    allpairs=[allpairs reshape(permute(pairTypes,[3,1,2]),...
        [size(corrTable.cellTypes{1},3),size(pairTypes,1)*size(pairTypes,2)])];
end
    
% first histogram our ratios
% we can recalculate our peaks here
newpeaks=0;
window=.1;

if newpeaks==1
    for i=1:size(pairIndex,2)
        [mypeak,myind]=max(pairIndex(abs(span)<window,i));
        newspan=span(abs(span)<window);
        %if mypeak>2
            
            allstats(1,i)=mypeak;
            allstats(2,i)=newspan(myind);
        %else
        %    allstats(1,i)=nan;
        %    allstats(2,i)=nan;
        %end
    end
end

figure;
%minus values means that the first vector (PFC) leads
% first get a waterfall
subplot(2,3,[1 4]);
okpairs=allstats(2,allstats(1,:)>2);
okgrams=pairIndex(:,allstats(1,:)>2);
[~,a]=sort(okpairs);
imagesc(span,1:length(okpairs),zscore(okgrams(:,a))');
hold on; plot([0 0],[1 length(a)],'k');
yyaxis right;
% mean
%plot(span,mean(okgrams')./max(mean(okgrams')),'w');

% or histogram
[cts,bins]=histcounts(okpairs,-window:.005:window);
plot(mean([bins(2:end); bins(1:end-1)]),cts./max(cts),'w-o');
[~,p2]=ttest(okpairs);
title(sprintf('All Cross Region \n p=%.2e p=%.2e',signrank(okpairs),p2));
xlabel('PFC leads <-> CA1 leads');

subinds=[2 3 5 6];
mytitles={'PFC INs <-> CA1 INs','PFC PYRs <-> CA1 INs',...
    'PFC INs <-> CA1 PYRs','PFC PYRs <-> CA1 PYRs'};

for i=1:4
% now for the four subtypes
subplot(2,3,subinds(i));
okpairs=allstats(2,allstats(1,:)>2 & allpairs==i-1);
okgrams=pairIndex(:,allstats(1,:)>2 & allpairs==i-1);
[~,a]=sort(okpairs);
imagesc(span,1:length(okpairs),zscore(okgrams(:,a))');
hold on; plot([0 0],[1 length(a)],'k');
yyaxis right;
%plot(span,mean(okgrams')./max(mean(okgrams')),'w');

[cts,bins]=histcounts(okpairs,-window:.005:window);
plot(mean([bins(2:end); bins(1:end-1)]),cts./max(cts),'w-o');
[~,p2]=ttest(okpairs);
title(sprintf('%s \n p=%.2e p=%.2e',mytitles{i},signrank(okpairs),p2));
xlabel('PFC leads <-> CA1 leads');
end

linkaxes(get(gcf,'Children'),'x');
xlim([-.15 .15]);
%%


figure;
histogram(allstats(3,:),50)
betapairs=pairIndex(:,allstats(3,:)<-.1);
betapeaks=allstats(2,allstats(3,:)<-.1);
[~,a]=sort(betapeaks);
title('theta-beta ratio');
subplot(1,3,2);
imagesc(span,1:length(a),zscore(betapairs(:,a))');
hold on; plot([0 0],[1 length(a)],'k');
yyaxis right;
[a,b]=histcounts(betapeaks,15);
%plot(mean([b(2:end); b(1:end-1)])',a,'r');
title(sprintf('p=%.2e',signrank(betapeaks)));


resppairs=pairIndex(:,allstats(3,:)>-.1);
resppeaks=allstats(2,allstats(3,:)>-.1);
[~,a]=sort(resppeaks);
subplot(1,3,3);
imagesc(span,1:length(a),zscore(resppairs(:,a))');
hold on; plot([0 0],[1 length(a)],'k');
[a,b]=histcounts(resppeaks,-.08:.01:.08);
yyaxis right;
%plot(mean([b(2:end); b(1:end-1)])',a,'r');
title(sprintf('p=%.2e',signrank(resppeaks)));
%%
% what about direction => frequency?

figure;
betapairs=pairIndex(:,allstats(2,:)<0);
betapeaks=allstats(2,allstats(2,:)<0);
betaratio=allstats(3,allstats(2,:)<0);
[~,a]=sort(betapeaks);
subplot(2,2,1);
imagesc(span,1:length(a),zscore(betapairs(:,a))');
hold on; plot([0 0],[1 length(a)],'k');
subplot(2,2,3);
histogram(betaratio,15);
%plot(mean([b(2:end); b(1:end-1)])',a,'r');
%title(sprintf('p=%.2e',signrank(betapeaks)));
subplot(2,2,2);

resppairs=pairIndex(:,allstats(2,:)>0);
resppeaks=allstats(2,allstats(2,:)>0);
respratio=allstats(3,allstats(2,:)>0);

[~,a]=sort(resppeaks);
imagesc(span,1:length(a),zscore(resppairs(:,a))');
hold on; plot([0 0],[1 length(a)],'k');
subplot(2,2,4);
histogram(respratio,15);

%plot(mean([b(2:end); b(1:end-1)])',a,'r');
title(sprintf('p=%.2e',ranksum(betaratio,respratio)));






%%
%cf=[.3 0 2 8 .1 .01 0];
hfig=figure;
subplot(2,1,1);
slow=(cf(1)*(cos(2*pi*cf(4)*(x+cf(7)))+1)+cf(2)).*exp(-abs((x+cf(7)))/cf(5));
fast=cf(3)*exp(-((x+cf(7)).^2)/cf(6)^2);

b=bar(span,normcurve,'FaceColor','flat','EdgeColor','flat');
b.CData=[.7 .7 .7]; hold on;
plot(x,slow,'k'); plot(x,fast,'r'); plot(x,slow+fast,'b');
legend('data','slowRR','fast','fit');

subplot(2,1,2);
slow=(cf2(1)*(cos(2*pi*cf2(4)*(x+cf2(7)))+1)+cf2(2)).*exp(-abs((x+cf2(7)))/cf2(5));
fast=cf2(3)*exp(-((x+cf2(7)).^2)/cf2(6)^2);
b=bar(span,normcurve,'FaceColor','flat','EdgeColor','flat');
b.CData=[.7 .7 .7]; hold on;
plot(x,slow,'k'); plot(x,fast,'r'); plot(x,slow+fast,'b');
legend('data','slowBeta','fast','fit');




        
title(sprintf('%.5f, %.5f, %.5f',cf(2), cf(3), cf(4)));
%%
%
%
%
%
%  Recapitulate claires cell cell analysis exactly
%
%
%
%
%%
%{
the key behind this analysis is counting the number of spikes at each time
binned delay.  So instead of xcorr being the engine, its nspikes that occur
at that time delay.  This is basically the same as using the other method.

The other key is that I jitter all spikes over a 50 msec bin so that the
overall crown of the cross correlogram is reflected in the null.  This
crown is destroyed if you shuffle trials (overall firing rate fluctuates by
trial, so mixing causes high fr trials to corr to low fr trials) and if you
circshift within trial, because the within-trial rate changes are destroyed
(either increasing overall, or decreasing overall)

Here i also use a very long odor duration, which could be shifting back
towards pre-sampling period.  I will try to shorten this duration and see
if the data still fall out, but i get the impression the duration is just
too short to do it this way


%}
%% Gather the cell-cell xcorrellograms, find crucial stats

% matt jones and matt wilson, 1005 plos theta rhythms coordinate hc pfc
% interactions in spatial memory task

% basically they measure the peak from +200msec to -200 msec and do a
% difference score

% I also want to measure the amount of beta or theta in the xcgram for c
% and incorrect

% use the get ACinfo to measure beta and theta in autocorr, do this for
% correct first just to see what most pairs do, then see if it's different
% for incorrect trials (will need to do some stats here)


% to my knowledge here are claires numbers:
% tbin=.01; corrspan=.2; skernel=2; % bins 

% Here I correct for edge effects.... by taking a trialwise correlation
tbin=.0025; corrspan=.3;
span=-corrspan:tbin:corrspan;
regions={'PFC','CA1'};
nboots=200;  plotIT=0;  odorDur=0.5;

pairIndex=[0 0; 1 0; 0 1; 1 1];

varNames={'xcorrgrams','statstable','celltypes','rawgrams'};

corrTable=table('Size',[length(SuperRat),4],'VariableNames',varNames,'VariableTypes',repmat({'cell'},1,4));
allclock=tic;
for i=1:length(SuperRat)
    sessclock=tic;
       
    % start with all cells and PFC is first region (as in regions)
    PFCunits=find(cellfun(@(a) contains(a,regions{1}), {SuperRat(i).units.area}) & ...
        cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
    
    CA1units=find(cellfun(@(a) contains(a,regions{2}), {SuperRat(i).units.area}) & ...
        cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
    
    % all trials, not c vs i yet
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    oktrials=diff(odorTimes,1,2)<2.5 & diff(odorTimes,1,2)>.5; % from 0.5 to 2 second samples
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10; % from 0.5 to 2 second samples
    
    % use only okay trials
    odorTimes=odorTimes(oktrials,:);
    trialCorrect=trialCorrect(oktrials,:);
    
    % cast all odorTimes to 0.5 seconds
    odorTimes(:,1)=odorTimes(:,2)-odorDur;
    Ctrials=odorTimes(trialCorrect==1,:);
    Itrials=odorTimes(trialCorrect==0,:);
    % corrmat will be big, it'll be a 3d, with i pfc cells by j HPC cells by k corrvals
    corrMatrix=nan(length(PFCunits),length(CA1units),length(span));
    rawMatrix=corrMatrix;
    corrStats=nan(length(PFCunits),length(CA1units),3); 
    
    % celltypetable- 0=IN IN; 1=PYR-IN; 2=IN-PYR; 3=PYR-PYR
    cellTypes=nan(length(PFCunits),length(CA1units));
    
    for j=1:length(PFCunits)
        % pspikes are the raw ts, pspk is how you bootstrap
        [Pspikes,~,~,~,~,Pspk]=event_spikes(SuperRat(i).units(PFCunits(j)).ts(:,1),...
            Ctrials(:,1),0,Ctrials(:,2)-Ctrials(:,1));
        Ptype=contains(SuperRat(i).units(PFCunits(j)).type,'pyr');
        
        for k=1:length(CA1units)
            [Cspikes,~,~,~,~,Cspk]=event_spikes(SuperRat(i).units(CA1units(k)).ts,...
                Ctrials(:,1),0,Ctrials(:,2)-Ctrials(:,1));
            Ctype=contains(SuperRat(i).units(CA1units(k)).type,'pyr');
            
            cellTypes(j,k)=find(sum([Ptype Ctype]==pairIndex,2)==2)-1;
            % count trials where both cells fire
            activetrials=cellfun(@(a,b) ~isempty(a) & ~isempty(b), Pspk, Cspk);
            
            if sum(activetrials)>=4
                
                trCorr=[];
                % this does all the spikes, lots of unnecessary data
                %spikediffs=Pspikes-Cspikes';
                spikediffs=cell2mat(cellfun(@(a,b) linearize(a-b'),Pspk,Cspk,'UniformOutput',false)'); 
                for td=1:length(span)
                    trCorr(td)=sum(abs(spikediffs(:)-span(td))<=tbin/2);
                end
                
                % for boot, remake the trialwise by adding starttimes to
                % Pspk
                bootCorr=nan(nboots,length(span));
                %mixup=(rand(length(Pspikes),nboots)-.5)*.1; % center around zero, span 100 msec
                %Pmix=mixup+Pspikes;
                
                for bt=1:nboots
                    % add or subtract up to 50 msec to each spike
                    Pmix=cellfun(@(a) a+(rand(length(a),1)-.5)*.1, Pspk,'UniformOutput',false);
                    spikediffs=cell2mat(cellfun(@(a,b) linearize(a-b'),Pmix,Cspk,'UniformOutput',false)'); 
                    %{
                    trCorrBt=[]; Pmix={};
                    for tr=1:length(Pspk)
                        %Pmix{tr}=mixup{tr}+Ctrials(tr,1);
                        jitter=rand(length(Pspk{tr}),1)*.05;
                        Pmix{tr}=Pspk{tr}+jitter+Ctrials(tr,1);
                    end
                    % again, if small negative, PFC leads
                   
                    spikediffs=cell2mat(Pmix')-Cspikes';
                    %}
                    for td=1:length(span)
                        bootCorr(bt,td)=sum(abs(spikediffs-span(td))<=tbin/2);
                    end
                end
                              
                mu=mean(bootCorr,1,'omitnan'); sig=std(bootCorr,1,'omitnan');
                %bootZ=(trCorr-mu)./sig;
                bootZ=mean(trCorr>bootCorr);
                %{  
                figure; subplot(2,1,1);
                plot(span,trCorr); hold on;                
                patch([span fliplr(span)], [mu+sig*2 fliplr(mu-sig*2)],[.7 .7 .7],'LineStyle','none','FaceAlpha',.3);
                subplot(2,1,2); plot(span,bootZ);   
                %}
                normcurve=rescale(trCorr);                
                corrMatrix(j,k,:)=bootZ;
                rawMatrix(j,k,:)=normcurve;
                [corrStats(j,k,1),peakind]=max(bootZ); corrStats(j,k,2)=span(peakind);
                %corrMatrix(j,k,:)=normcurve;                
            end
        end
    end
    
    corrTable.xcorrgrams(i)={corrMatrix};
    corrTable.corrStats(i)={corrStats};
    corrTable.cellTypes(i)={cellTypes};
    corrTable.rawgrams(i)={rawMatrix};
    % on first blush, it looks like specific PFC units will assume
    % peaks in their xcorr with many CA1 units during odor sampling
    % she will plot the mean of all these
    fprintf('Sess %d, rat %s %d done in %d, prob %d more\n',i,SuperRat(i).name,...
    SuperRat(i).daynum,round(toc(sessclock)), round(toc(allclock)/i*(length(SuperRat)-i)));
end

% corrStats will be empty
%%

corrTable(cellfun(@(a) isempty(a), corrTable.xcorrgrams),:)=[];
allGrams=[]; allstats=[]; allpairs=[];
for i=1:height(corrTable)
    %thisCorr=corrTable.xcorrgrams{i};
    thisCorr=corrTable.rawgrams{i};
    thisStats=corrTable.corrStats{i};
    pairTypes=corrTable.cellTypes{i};
    allstats=[allstats reshape(permute(thisStats,[3,1,2]),...
        [3,size(thisStats,1)*size(thisStats,2)])];
    allGrams=[allGrams reshape(permute(thisCorr,[3,1,2]),...
        [size(corrTable.xcorrgrams{1},3),size(thisCorr,1)*size(thisCorr,2)])];
    allpairs=[allpairs reshape(permute(pairTypes,[3,1,2]),...
        [size(corrTable.cellTypes{1},3),size(pairTypes,1)*size(pairTypes,2)])];
end
    
%allGrams(:,any(isnan(allGrams),1))=[];
% first histogram our ratios
% we can recalculate our peaks here
newpeaks=1;
window=.1;
winbin=.005;
skernel=1;
%allGrams=zscore(SmoothMat2(allGrams,[0 8],skernel));
windowedges=-window-winbin/2:winbin:window+winbin/2;
if newpeaks==1
    for i=1:size(allGrams,2)
        [mypeak,myind]=max(allGrams(abs(span)<window,i));
        newspan=span(abs(span)<window);

        allstats(1,i)=mypeak;
        allstats(2,i)=newspan(myind);
        allstats(3,i)=myind;
        if any(isnan(allGrams(:,i))) % ||mypeak<2  
            allstats(1,i)=nan;
            allstats(2,i)=nan;
            allstats(3,i)=nan;
        end

            
    end
end
allGrams=SmoothMat2(allGrams,[0 8],skernel);
figure;
%minus values means that the first vector (PFC) leads
% first get a waterfall
corrthresh=.99;
subplot(2,3,[1 4]);
okpairs=allstats(2,allstats(1,:)>corrthresh);
okgrams=allGrams(:,allstats(1,:)>corrthresh);
[~,a]=sort(okpairs);
if Zfirst==1
    imagesc(span,1:length(okpairs),okgrams(:,a)');
else
    imagesc(span,1:length(okpairs),zscore(okgrams(:,a))'); % smooth later or
end

hold on; plot([0 0],[1 length(a)],'k');
yyaxis right;
% mean
%plot(span,mean(okgrams')./max(mean(okgrams')),'w');

% or histogram
[cts,bins]=histcounts(okpairs,windowedges); % no smooth here
wincts=okpairs'-(windowedges(2:end)-winbin/2);
cts=sum(abs(wincts)<=winbin*2);
plot(mean([bins(2:end); bins(1:end-1)]),cts./max(cts),'k-*',...
    'MarkerSize',4,'LineWidth',1.2);
[~,p2]=ttest(okpairs);
title(sprintf('All Cross Region \n x=%.3f p=%.2e p=%.2e',...
    mean(okpairs),signrank(okpairs(abs(okpairs)<.05)),p2));
xlabel('PFC leads <-> CA1 leads');
ylim([0 4]); set(gca,'YTick',[0 .5 1]);

subinds=[2 3 5 6];
mytitles={'PFC INs <-> CA1 INs','PFC PYRs <-> CA1 INs',...
    'PFC INs <-> CA1 PYRs','PFC PYRs <-> CA1 PYRs'};

for i=1:4
% now for the four subtypes
subplot(2,3,subinds(i));
okpairs=allstats(2,allstats(1,:)>corrthresh & allpairs==i-1);
okgrams=allGrams(:,allstats(1,:)>corrthresh & allpairs==i-1);
[~,a]=sort(okpairs);
if Zfirst==1
    imagesc(span,1:length(okpairs),okgrams(:,a)');
else
    imagesc(span,1:length(okpairs),zscore(okgrams(:,a))'); % smooth later or
end
hold on; plot([0 0],[1 length(a)],'k');
yyaxis right;
%plot(span,mean(okgrams')./max(mean(okgrams')),'w');

[cts,bins]=histcounts(okpairs,windowedges);
wincts=okpairs'-(windowedges(2:end)-winbin/2);
cts=sum(abs(wincts)<=winbin*2);
%plot(mean([bins(2:end); bins(1:end-1)]),SmoothMat2(cts./max(cts),[1 5],1),'k-*',...
%    'MarkerSize',4,'LineWidth',1.2);
plot(mean([bins(2:end); bins(1:end-1)]),cts./max(cts),'k-*',...
    'MarkerSize',4,'LineWidth',1.2);
[~,p2]=ttest(okpairs);
title(sprintf('%s \n x=%.3f n=%d p=%.2e',mytitles{i},mean(okpairs),length(okpairs),signrank(okpairs)));
xlabel('PFC leads <-> CA1 leads');
ylim([0 4]); set(gca,'YTick',[0 .5 1]);
end

linkaxes(get(gcf,'Children'),'x');
xlim([-window window]);
sgtitle(sprintf('Bins=%f, 1 bin kernel',tbin));

%% align all to center

figure;
%minus values means that the first vector (PFC) leads
% first get a waterfall
corrthresh=.7;
subplot(2,3,[1 4]);
okpairs=allstats(2,allstats(1,:)>corrthresh);
okgrams=allGrams(:,allstats(1,:)>corrthresh);
okoffsets=allstats(3,allstats(1,:)>corrthresh);

% first circshift these
shiftnum=okoffsets-ceil(length(newspan)/2);
for i=1:size(okgrams,2)
    okgrams2(:,i)=circshift(okgrams(:,i),-shiftnum(i));
end

[~,a]=sort(okpairs);

%imagesc(span,1:length(okpairs),zscore(okgrams(:,a))');
imagesc(span,1:length(okpairs),okgrams2(:,a)');

hold on; plot([0 0],[1 length(a)],'k');
yyaxis right;
% mean
plot(span,mean(okgrams2','omitnan')./max(mean(okgrams2','omitnan')),'k');
ylim([-1 4])

for i=1:4
% now for the four subtypes
subplot(2,3,subinds(i));
okpairs=allstats(2,allstats(1,:)>corrthresh & allpairs==i-1);
okgrams=allGrams(:,allstats(1,:)>corrthresh & allpairs==i-1);
okoffsets=allstats(3,allstats(1,:)>corrthresh & allpairs==i-1);

% first circshift these
shiftnum=okoffsets-ceil(length(newspan)/2);
for k=1:size(okgrams,2)
    okgrams2(:,k)=circshift(okgrams(:,k),-shiftnum(k));
end

[~,a]=sort(okpairs);
imagesc(span,1:length(okpairs),zscore(okgrams2(:,a))');
%imagesc(span,1:length(okpairs),okgrams(:,a)');
hold on; plot([0 0],[1 length(a)],'k');
yyaxis right;
%plot(span,mean(okgrams')./max(mean(okgrams')),'w');

[cts,bins]=histcounts(okpairs,windowedges);
%  check smoothing kernel here
plot(mean([bins(2:end); bins(1:end-1)]),SmoothMat2(cts./max(cts),[1 5],1),'k-*',...
    'MarkerSize',4,'LineWidth',1.2);
[~,p2]=ttest(okpairs);
title(sprintf('%s \n x=%.3f p=%.2e p=%.2e',mytitles{i},mean(okpairs),signrank(okpairs),p2));
xlabel('PFC leads <-> CA1 leads');
ylim([-1 4]); set(gca,'YTick',[0 .5 1]);
end

linkaxes(get(gcf,'Children'),'x');
xlim([-.15 .15]);

%%
%
%
%
%  run analysis agnostic to region, parse regions later
%
%
%
%
%% Gather the cell-cell xcorrellograms, find crucial stats

% matt jones and matt wilson, 1005 plos theta rhythms coordinate hc pfc
% interactions in spatial memory task

% basically they measure the peak from +200msec to -200 msec and do a
% difference score

% I also want to measure the amount of beta or theta in the xcgram for c
% and incorrect

% use the get ACinfo to measure beta and theta in autocorr, do this for
% correct first just to see what most pairs do, then see if it's different
% for incorrect trials (will need to do some stats here)


% Here I correct for edge effects.... by taking a trialwise correlation
tbin=.0025;  gauss=.015;  corrspan=.3;
span=-corrspan:tbin:corrspan;
regions={'PFC','CA1'};
nboots=100;  plotIT=0;  odorDur=0.5;

% so it will be
% cellTypes- 0=IN IN; 1=PYR-IN; 2=IN-PYR; 3=PYR-PYR
% regionInds
pairIndex=[0 0; 1 0; 0 1; 1 1];

varNames={'odorgrams','rungrams','cellTypes','regionPairs'};
corrTable=table('Size',[length(SuperRat),4],'VariableNames',varNames,'VariableTypes',repmat({'cell'},1,4));
allclock=tic;
for i=1:length(SuperRat)
    sessclock=tic;
     
    % all trials, not c vs i yet
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    try
        % add odor dur to establish a buffer time (same as bayesian
        % analysis)
        runTimes=[SuperRat(i).trialdata.sniffend+odorDur SuperRat(i).trialdata.rewardstart];
    catch
        continue
    end
    oktrials=diff(odorTimes,1,2)<2.5 & diff(odorTimes,1,2)>.5; % from 0.5 to 2 second samples
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10; % from 0.5 to 2 second samples
    
    % use only okay trials
    odorTimes=odorTimes(oktrials,:);
    trialCorrect=trialCorrect(oktrials,:);
    runTimes=runTimes(oktrials,:);
    
    % cast all odorTimes to 0.5 seconds
    odorTimes(:,1)=odorTimes(:,2)-odorDur;
    Ctrials=odorTimes(trialCorrect==1,:); Cruns=runTimes(trialCorrect==1,:);
    Itrials=odorTimes(trialCorrect==0,:); Iruns=runTimes(trialCorrect==0,:);

    % corrmat will be big, it'll be a 3d, with i pfc cells by j HPC cells by k corrvals
    odorMatrix=nan(length(SuperRat(i).units),length(SuperRat(i).units),length(span)); % -25:25 bins
    odorStats=nan(length(SuperRat(i).units),length(SuperRat(i).units),3); % this will be center in msec, and theta-gamma ratio
    runMatrix=nan(length(SuperRat(i).units),length(SuperRat(i).units),length(span)); % -25:25 bins
    runStats=nan(length(SuperRat(i).units),length(SuperRat(i).units),3); % this will be center in msec, and theta-gamma ratio
    
    % cellTypes- 0=IN IN; 1=PYR-IN; 2=IN-PYR; 3=PYR-PYR
    % regionPairs 0=PFC-PFC; 1=CA1-PFC; 2=PFC-CA1; 3=CA1-CA1
    cellTypes=nan(length(SuperRat(i).units),length(SuperRat(i).units));
    regionPairs=cellTypes;
    
    
    
    for j=1:length(SuperRat(i).units)-1
        % pspikes are the raw ts, pspk is how you bootstrap
        [Pspikes,~,~,~,~,Pspk]=event_spikes(SuperRat(i).units(j).ts(:,1),...
            Ctrials(:,1),0,Ctrials(:,2)-Ctrials(:,1));
        [Prspikes,~,~,~,~,Prspk]=event_spikes(SuperRat(i).units(j).ts(:,1),...
            Cruns(:,1),0,Cruns(:,2)-Cruns(:,1));
        
        Ptype=contains(SuperRat(i).units(j).type,'pyr');
        Preg=contains(SuperRat(i).units(j).area,'CA1');
        
        for k=j+1:length(SuperRat(i).units)
            
            if SuperRat(i).units(j).tet~=SuperRat(i).units(k).tet
                [Cspikes,~,~,~,~,Cspk]=event_spikes(SuperRat(i).units(k).ts,...
                    Ctrials(:,1),0,Ctrials(:,2)-Ctrials(:,1));
                [Crspikes,~,~,~,~,Crspk]=event_spikes(SuperRat(i).units(k).ts,...
                    Cruns(:,1),0,Cruns(:,2)-Cruns(:,1));
                
                Ctype=contains(SuperRat(i).units(k).type,'pyr');
                Creg=contains(SuperRat(i).units(k).area,'CA1');
                
                cellTypes(j,k)=find(sum([Ptype Ctype]==pairIndex,2)==2)-1;
                regionTypes(j,k)=find(sum([Preg Creg]==pairIndex,2)==2)-1;
                % count trials where both cells fire
                activetrials=cellfun(@(a,b) ~isempty(a) & ~isempty(b), Pspk, Cspk);
                % for odor period
                if sum(activetrials)>=4
                    
                    trCorr=[]; spikediffs=Pspikes-Cspikes';
                    for td=1:length(span)
                        trCorr(td)=sum(spikediffs(:)>span(td)-tbin/2 & spikediffs(:)<=span(td)+tbin/2);
                    end
                    % for boot, remake the trialwise by adding starttimes to Pspk
                    bootCorr=[]; mixup=rand(length(Pspikes),nboots)*.05;
                    for bt=1:nboots
                        % add or subtract up to 50 msec to each spike
                        Pmix=Pspikes+mixup(:,bt);
                        spikediffs=Pmix'-Cspikes;
                        for td=1:length(span)
                            bootCorr(bt,td)=sum(spikediffs(:)>span(td)-tbin/2 & spikediffs(:)<=span(td)+tbin/2);
                        end
                    end
                    
                    mu=mean(bootCorr,1,'omitnan'); sig=std(bootCorr,1,'omitnan');
                    bootZ=(trCorr-mu)./sig;
                    odorMatrix(j,k,:)=bootZ;
                    %corrMatrix(j,k,:)=normcurve;
                    
                figure; subplot(2,2,1);
                plot(span,trCorr); hold on;
                patch([span fliplr(span)], [mu+sig*2 fliplr(mu-sig*2)],[.7 .7 .7],'LineStyle','none','FaceAlpha',.3);
                subplot(2,2,2); plot(span,bootZ);
                    
                    
                end
                % for runs
                activetrials=cellfun(@(a,b) ~isempty(a) & ~isempty(b), Prspk, Crspk);
                if sum(activetrials)>=4
                    
                    trCorr=[]; spikediffs=Prspikes-Crspikes';
                    for td=1:length(span)
                        trCorr(td)=sum(spikediffs(:)>span(td)-tbin/2 & spikediffs(:)<=span(td)+tbin/2);
                    end
                    % for boot, remake the trialwise by adding starttimes to Pspk
                    bootCorr=[]; mixup=rand(length(Prspikes),nboots)*.05;
                    for bt=1:nboots
                        % add or subtract up to 50 msec to each spike
                        Pmix=Prspikes+mixup(:,bt);
                        spikediffs=Pmix'-Crspikes;
                        for td=1:length(span)
                            bootCorr(bt,td)=sum(spikediffs(:)>span(td)-tbin/2 & spikediffs(:)<=span(td)+tbin/2);
                        end
                    end
                    
                    mu=mean(bootCorr,1,'omitnan'); sig=std(bootCorr,1,'omitnan');
                    bootZ=(trCorr-mu)./sig;
                    runMatrix(j,k,:)=bootZ;
                    %corrMatrix(j,k,:)=normcurve;
                    
                subplot(2,2,3);
                plot(span,trCorr); hold on;
                patch([span fliplr(span)], [mu+sig*2 fliplr(mu-sig*2)],[.7 .7 .7],'LineStyle','none','FaceAlpha',.3);
                subplot(2,2,4); plot(span,bootZ);

                    
                end
            end
        end
    end
    
    corrTable.odorGrams(i)={odorMatrix};
    corrTable.runGrams(i)={runMatrix};
    corrTable.cellTypes(i)={cellTypes};
    corrTable.regionPairs(i)={regionPairs};
    % on first blush, it looks like specific PFC units will assume
    % peaks in their xcorr with many CA1 units during odor sampling
    % she will plot the mean of all these
    fprintf('Sess %d, rat %s %d done in %d, prob %d more\n',i,SuperRat(i).name,...
    SuperRat(i).daynum,round(toc(sessclock)), round(toc(allclock)/i*(length(SuperRat)-i)));
end