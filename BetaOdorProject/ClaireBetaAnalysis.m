% ClaireBetaAnalysis

%% Comments block
%{


Claires algorithm:
1. find tet with most cells
2. get that tets beta filtered lfp
3. interp eeg times

4. get spikes in odor delivery period
        or epoch spikes based on when beta envelope is >2sd above mean,
        using a fwhm approach.
5. get phase of those spikes
6. calc mean phase, mean vector length, p value
7. bootstrap that dataset (could probably pull beta from different trials)

8. aggregate for pyrams and INs and pfc and hpc



%%%%% overall results %%%%
1. there are hpc and pfc pyrs and ins that cohere to rr and beta
2. the big one is that HPC cells cohere to beta
3. and these cells cohere less during error trials

Panels.
a. example beta, RR coherent neurons
b. nn cells in pfc/hpc cohering to  beta, rr same cross
    only within region coherence, stacked pyr in, x beta, rr, and color
    pair is ca1 pfc.   darker for INs

c.


%}

%clearvars -except SuperRat

% set region parameters here!!
regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')]; % BETA IS BLUE
types={'pyr','in'};
rhythm={'beta','resp'};
%% Overall beta coherence using whole odor period and just the correct trials
% we should be using stats only on local beta, but we can do all


plotIT=0;
    
for i=1:length(SuperRat)
    myclock=tic;
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorr=SuperRat(i).trialdata.CorrIncorr10;
    oktrials=trialCorr & diff(snifftimes,1,2)>=.5 & diff(snifftimes,1,2)<2.5;


    for k=1:length(regions)
        lfpSource=matfile(SuperRat(i).([regions{k} 'eegFile']));
        for rh=1:2
            myLFP=lfpSource.([rhythm{rh} 'continuous']);
        for j=1:length(SuperRat(i).units)
            myspikes=EpochCoords(SuperRat(i).units(j).ts(:,1),snifftimes(oktrials,:));
            % run on local beta, thats gonna be the best ref

            % insert here case if no viable LFP
            
            SuperRat(i).units(j).([rhythm{rh} 'mean'])(k)=nan;
            SuperRat(i).units(j).([rhythm{rh} 'MVL'])(k)=nan;
            SuperRat(i).units(j).([rhythm{rh} 'Ostat'])(k)=nan;
            SuperRat(i).units(j).([rhythm{rh} 'Rstat'])(k)=nan;

            if length(myspikes)>5 && ~isempty(myLFP)

                % if there are enough spikes to analyze, get the mean
                % phase, the mvl, the r test and O test
                phasesCorr=interp1(myLFP(:,1),myLFP(:,3),myspikes,'nearest');
                SuperRat(i).units(j).([rhythm{rh} 'mean'])(k)=circ_mean(phasesCorr);
                SuperRat(i).units(j).([rhythm{rh} 'MVL'])(k)=circ_r(phasesCorr);
                SuperRat(i).units(j).([rhythm{rh} 'Ostat'])(k)=circ_otest(phasesCorr);
                SuperRat(i).units(j).([rhythm{rh} 'Rstat'])(k)=circ_rtest(phasesCorr);

                % plot it out???
                if SuperRat(i).units(j).betaRstat(k)<.01 && plotIT
                    figure; ha=histogram([phasesCorr; phasesCorr+pi*2],linspace(-pi,pi*3,10+ceil(sqrt(length(phasesCorr)))),...
                        'Normalization','probability');
                    set(ha,'LineStyle','none','FaceColor',colors(k,:));
                    set(gca,'XTick',[pi.*[-1:4]],'XTickLabel',{'-2\pi','','0','','2\pi'});
                    ylabel('Probability of Spike'); xlabel('Beta Phase'); box off;
                    title(sprintf(' %s cell in %s, %s lfp %s \n n=%.f, mean %.2f, mvl=%.2f, p=%.2e',...
                        SuperRat(i).units(j).type, SuperRat(i).units(j).area, rhythm{rh}, regions{k}, length(phasesCorr),...
                        SuperRat(i).units(j).betamean(k),SuperRat(i).units(j).betaMVL(k),SuperRat(i).units(j).betaRstat(k)));
                    uiwait;
                end
            end
        end
    end
    fprintf('Session %d, rat %s day %d done in %d seconds \n',i, SuperRat(i).name,SuperRat(i).daynum, round(toc(myclock)));
end

end




%% Compare strength of modulation for correct and incorrect trials

% she did rayleigh z, but that is a shitty metric.  I would like to measure
% the mean resultant vector length.  For that, I'll bootstrap or jackknife
% a mean and then measure the difference

%{
200 boots, and I basically alwatys downsample the number of spikes in
correct trials to the number in the incorrect trials.  For the larger
datset though, I bootstrap a mean coherence from 200 samples.

I also threshold for at least 10 spikes in the lesser of the two conditions
Spanning this parameter does little to change the outcome

and finally i interpolate to the nearest phase timestamp so that it doesnt
inproperly wrap the phases if its around the bend

the betadprime isnt in the figure, this is just to convince myself that on
a single cell basis, these coherence values actually meaningfully differ


%}

%% this is beta first
nboots=200; bootclock=tic;
rstream = RandStream('dsfmt19937','Seed',16);
RandStream.setGlobalStream(rstream);

for i=1:length(SuperRat)
    myclock=tic;
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    oktimes=diff(snifftimes,1,2)>=.5 & diff(snifftimes,1,2)<2.5;
    snifftimes=snifftimes(oktimes,:);
    trialCorr=SuperRat(i).trialdata.CorrIncorr10(oktimes);


    h=waitbar(0,sprintf('Starting beta animal %s %d', SuperRat(i).name,SuperRat(i).daynum));
    for j=1:length(SuperRat(i).units)
        [~,~,~,~,myspikes]=event_spikes(SuperRat(i).units(j).ts(:,1),snifftimes(:,1),0,snifftimes(:,2)-snifftimes(:,1));
        
        for k=1:length(regions)
            lfpSource=matfile(SuperRat(i).([regions{k} 'eegFile']));
            myLFP=lfpSource.(['betacontinuous']);
            % insert here case if no viable LFP
            Cspikes=cell2mat(myspikes(trialCorr==1)');
            Ispikes=cell2mat(myspikes(trialCorr~=1)');
            % pull nspikes for the smaller of the two
            minspkct=round(min([length(Cspikes) length(Ispikes)]));
            
            if minspkct>10 && ~isempty(myLFP)
                Cphases=interp1(myLFP(:,1),myLFP(:,3),Cspikes,'nearest');
                Iphases=interp1(myLFP(:,1),myLFP(:,3),Ispikes,'nearest');

                % there are two thinning procedures I can use- one- thin spikes
                % regardless of trial, or 2 thin trials and therefore spikes
                
                % if there are enough spikes to analyze
                bootC=nan(1,nboots); % allocate the two distributions
                bootI=nan(1,nboots);
                parfor boot=1:nboots
                    bootC(boot)=circ_r(datasample(Cphases,minspkct,'Replace',false));
                    bootI(boot)=circ_r(datasample(Iphases,minspkct,'Replace',false));
                end
                % get a z value of the single one vs the multiple
                    SuperRat(i).units(j).betaPhaseDprime(k,:)=[mean(bootC,'omitnan') mean(bootI,'omitnan')...
                        (mean(bootC)-mean(bootI))/std(unique([bootC bootI]))];
                
            else
                SuperRat(i).units(j).betaPhaseDprime(k,:)=[nan nan nan];
            end
            
            % nan out our variables here
        end
        waitbar(j/length(SuperRat(i).units),h,sprintf('Running Unit %d of %d',j,length(SuperRat(i).units)));
    end
    close(h);
    fprintf('%s day %d took %d minutes, probably %d minutes left \n', SuperRat(i).name,...
        SuperRat(i).daynum, round(toc(bootclock)/i/60), round((toc(bootclock)/i/60)*(length(SuperRat)-i)))
end

%% this is for rr
rstream = RandStream('dsfmt19937','Seed',16);
RandStream.setGlobalStream(rstream);
nboots=200; bootclock=tic;
for i=13:length(SuperRat)
    myclock=tic;
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    oktimes=diff(snifftimes,1,2)>=.5 & diff(snifftimes,1,2)<2.5;
    snifftimes=snifftimes(oktimes,:);
    trialCorr=SuperRat(i).trialdata.CorrIncorr10(oktimes);

    h=waitbar(0,sprintf('Starting animal %s %d', SuperRat(i).name,SuperRat(i).daynum));
    for j=1:length(SuperRat(i).units)
        [~,~,~,~,myspikes]=event_spikes(SuperRat(i).units(j).ts(:,1),snifftimes(:,1),0,snifftimes(:,2)-snifftimes(:,1));
        
        for k=1:length(regions)
            lfpSource=matfile(SuperRat(i).([regions{k} 'eegFile']));
            myLFP=lfpSource.(['respcontinuous']);            % insert here case if no viable LFP
            Cspikes=cell2mat(myspikes(trialCorr==1)');
            Ispikes=cell2mat(myspikes(trialCorr~=1)');
            % pull nspikes for the smaller of the two
            minspkct=round(min([length(Cspikes) length(Ispikes)]));
            
            if minspkct>10 && ~isempty(myLFP)
                Cphases=interp1(myLFP(:,1),myLFP(:,3),Cspikes,'nearest');
                Iphases=interp1(myLFP(:,1),myLFP(:,3),Ispikes,'nearest');

                % there are two thinning procedures I can use- one- thin spikes
                % regardless of trial, or 2 thin trials and therefore spikes
                
                % if there are enough spikes to analyze
                bootC=nan(1,nboots); % allocate the two distributions
                bootI=nan(1,nboots);
                parfor boot=1:nboots
                    bootC(boot)=circ_r(datasample(Cphases,minspkct,'Replace',false));
                    bootI(boot)=circ_r(datasample(Iphases,minspkct,'Replace',false));
                end
                % get a z value of the single one vs the multiple
                
                    SuperRat(i).units(j).respPhaseDprime(k,:)=[mean(bootC,'omitnan') mean(bootI,'omitnan')...
                        (mean(bootC)-mean(bootI))/std(unique([bootI bootC]))]; % workaround for repeats!
            else
                SuperRat(i).units(j).respPhaseDprime(k,:)=[nan nan nan];
            end
            
            % nan out our variables here
        end
        waitbar(j/length(SuperRat(i).units),h,sprintf('Running Unit %d of %d',j,length(SuperRat(i).units)));
    end
    close(h);
    fprintf('%s day %d took %d minutes, probably %d minutes left \n', SuperRat(i).name,...
        SuperRat(i).daynum, round(toc(bootclock)/i/60), round((toc(bootclock)/i/60)*(length(SuperRat)-i)))
end

save('E:\Brandeis datasets\Claire Data\ClaireData-2022-07-21.mat','SuperRat','-v7.3')

%%
%
%
%
%
%        Plot functions for both beta and resp
%
%
%
%
%

% some indexing on cell responsivity and selectivity
%{
 regcoders=allcells(strcmpi({allcells.area},regions{i}) &...
            cellfun(@(a) a{1,2}>pcrit,{allcells.OdorSelective}) &...
            cellfun(@(a) a(3)<.05,{allcells.taskResponsive}) & ...
            strcmpi({allcells.type},type{t}));
%}

% need to run first block of claireRR analysis here
%% how many coherent cells in each region?
% now replicate claires analyses:
% local coherence and remote (cell in region 1, lfp region2)

% first, she shows a few beta locked cells....
% all odor responsive units in a given region

% for this analysis I think i'm just going to use sessions with all three
% regions present

pcrit=.05/3; % bonferonni correct, each cell is used 3 times

for celltype=1:2
clear respcount betacount;
for rh=1:2
    regions={'PFC','CA1'};
    regions2={'PFC','CA1','OB'};
    varTypes=repmat({'double'},1,6);
    allcount=table('size',[3 6],'VariableTypes',varTypes,'VariableNames',{'PFCtot','PFClocked','PFCpct','CA1tot','CA1locked','CA1pct'},...
        'rownames',{'PFCLFP','CA1LFP','OBLFP'});
    %contains({SuperRat(i).units.type},'in') & ...
    
    for i=1:length(SuperRat)
        for j=1:length(regions) % cell region
            mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
                cellfun(@(a) a(3)<.05, {SuperRat(i).units.taskResponsive}) &...
                contains({SuperRat(i).units.type},types{celltype}));
            for k=1:length(regions2) % lfp region
                allcount.([regions{j} 'tot'])(k)=allcount.([regions{j} 'tot'])(k)+length(mycells);
                % add the specific lfp rstat you want here (a(k))
                allcount.([regions{j} 'locked'])(k)=allcount.([regions{j} 'locked'])(k)...
                    +sum(cellfun(@(a) a(k),{mycells.([rhythm{rh} 'Rstat'])})<pcrit);
            end
        end
    end
    allcount.PFCpct=allcount.PFClocked./allcount.PFCtot;
    allcount.CA1pct=allcount.CA1locked./allcount.CA1tot;
    
    eval([types{celltype} rhythm{rh} 'coherence=allcount;']);
    eval(['openvar ' types{celltype} rhythm{rh} 'coherence']);
end
end

%% now for the nooverlap stacked bar
% i want to know each cells type, its origin, and then the beta and rr
% pvalues
pcrit=0.05/3;
allCells=cell2mat({SuperRat.units});
okcells=cell2mat(cellfun(@(a) a, {allCells.taskResponsive}, 'UniformOutput', false)');
allCells=allCells(okcells(:,3)<.05); % only take task responsive cells
cellTable=array2table(contains({allCells.area},'CA1')',"VariableNames",{'isCA1'});
cellTable.isPFC=contains({allCells.area},'PFC')';
cellTable.isPYR=contains({allCells.type},'pyr')';
cellTable.betaCoh=cellfun(@(a) any(a<pcrit), {allCells.betaRstat})';
[~,temp]=min(cell2mat(cellfun(@(a) a', {allCells.betaRstat}, 'UniformOutput', false)));
cellTable.betaPref=temp';
cellTable.respCoh=cellfun(@(a) any(a<pcrit), {allCells.respRstat})';
[~,temp]=min(cell2mat(cellfun(@(a) a', {allCells.respRstat}, 'UniformOutput', false)));
cellTable.respPref=temp';
cellTable.odorP=cellfun(@(a) a.issig(1), {allCells.OdorSelective})';
cellTable.odorSI=cellfun(@(a) a.score(1), {allCells.OdorSelective})';
% now make stacked bar on this:

bary=[accumarray(cellTable.betaPref(cellTable.betaCoh==1 & cellTable.isCA1==1 & cellTable.isPYR==1),1,[3, 1],@sum)'/sum(cellTable.isCA1==1 & cellTable.isPYR==1);...
   accumarray(cellTable.betaPref(cellTable.betaCoh==1 & cellTable.isCA1==1 & cellTable.isPYR==0),1,[3, 1],@sum)'/sum(cellTable.isCA1==1 & cellTable.isPYR==0);...
   accumarray(cellTable.betaPref(cellTable.betaCoh==1 & cellTable.isPFC==1 & cellTable.isPYR==1),1,[3, 1],@sum)'/sum(cellTable.isPFC==1 & cellTable.isPYR==1);...
   accumarray(cellTable.betaPref(cellTable.betaCoh==1 & cellTable.isPFC==1 & cellTable.isPYR==0),1,[3, 1],@sum)'/sum(cellTable.isPFC==1 & cellTable.isPYR==0)];

   figure;
   subplot(1,2,1);
barx=[.9 1.1 1.4 1.6];
bar(barx,bary,'stacked');
set(gca,'XTick', [.9 1.1 1.4 1.6],'XTickLabel',{'CA1 Pyr','CA1 IN','PFC Pyr','PFC IN'});

hold on; plot([0.7 1.9],[pcrit pcrit],'k--')
legend({'PFC beta','CA1 beta','OB beta'})
title('Beta Coherence')

bary2=[accumarray(cellTable.respPref(cellTable.respCoh==1 & cellTable.isCA1==1 & cellTable.isPYR==1),1,[3, 1],@sum)'/sum(cellTable.isCA1==1 & cellTable.isPYR==1);...
   accumarray(cellTable.respPref(cellTable.respCoh==1 & cellTable.isCA1==1 & cellTable.isPYR==0),1,[3, 1],@sum)'/sum(cellTable.isCA1==1 & cellTable.isPYR==0);...
   accumarray(cellTable.respPref(cellTable.respCoh==1 & cellTable.isPFC==1 & cellTable.isPYR==1),1,[3, 1],@sum)'/sum(cellTable.isPFC==1 & cellTable.isPYR==1);...
   accumarray(cellTable.respPref(cellTable.respCoh==1 & cellTable.isPFC==1 & cellTable.isPYR==0),1,[3, 1],@sum)'/sum(cellTable.isPFC==1 & cellTable.isPYR==0)];

subplot(1,2,2);
barx=[.9 1.1 1.4 1.6];
bar(barx,bary2,'stacked');
set(gca,'XTick', [.9 1.1 1.4 1.6],'XTickLabel',{'CA1 Pyr','CA1 IN','PFC Pyr','PFC IN'});

hold on; plot([0.7 1.9],[pcrit pcrit],'k--')
legend({'PFC RR','CA1 RR','OB RR'})
title('Respiratory Coherence')
linkaxes(get(gcf,'Children'));

% and now to tabulate the totals:

%% figure S5b, c, d, and e
% b: beta rr and overlap vs chance for CA1 and PFC PYRS
figure;
for i=1:2
bary=[mean(cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i));...
    mean(cellTable.respCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i));...
    mean(cellTable.respCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i) &...
    cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i)) ;...
    mean(cellTable.respCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i))*...
    mean(cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i));...

    mean(cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==2-i));...
    mean(cellTable.respCoh(cellTable.isPFC==1 & cellTable.isPYR==2-i));...
    mean(cellTable.respCoh(cellTable.isPFC==1 & cellTable.isPYR==2-i) &...
    cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==2-i)) ;...
    mean(cellTable.respCoh(cellTable.isPFC==1 & cellTable.isPYR==2-i))*...
    mean(cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==2-i))];


subplot(2,3,i);
bar([.9 1.1 1.4 1.6 2.4 2.6 2.9 3.1],bary');
set(gca,'XTick',[]), xlabel('CA1                    PFC');
ylabel(sprintf('Percentage of Task- \nResponsive %s Cells',types{i}))
% need to errorbar the grey bars (4 and 8)
% these can come from boostraps or just the bernoulli 95% likelihoods
% so binofit is the 95% conf, and binocdf is the likelihood
[~,pci]=binofit(round(bary(4)*sum(cellTable.isCA1==1 & cellTable.isPYR==1)),sum(cellTable.isCA1==1 & cellTable.isPYR==1));
[~,pci(2,:)]=binofit(round(bary(8)*sum(cellTable.isPFC==1 & cellTable.isPYR==1)),sum(cellTable.isPFC==1 & cellTable.isPYR==1));
pci=abs(pci-bary([4 8]));
hold on; errorbar([1.6 3.1],bary([4 8]),pci(:,1),pci(:,2),'k.','CapSize',0)
title('Beta and RR overlap');
% are the conditional probabilities the same too?
mean(cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i & cellTable.respCoh))
mean(cellTable.respCoh(cellTable.isCA1==1 & cellTable.isPYR==2-i & cellTable.betaCoh))
end
subplot(2,3,3)
bary=[mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==1));...
    mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==0));...
    mean(cellTable.odorP(cellTable.isCA1==0 & cellTable.isPYR==1));...
    mean(cellTable.odorP(cellTable.isCA1==0 & cellTable.isPYR==0))];
bar([.9 1.1 1.4 1.6],bary');
set(gca,'XTick',[]), xlabel(sprintf('Pyr IN         Pyr IN \nCA1                    PFC'));
ylabel(sprintf('Percentage of Task- \nResponsive Pyramidal Cells'));
ylim([0 .5])


% c: beta locked vs choice selective pyramidal cells
bary=[mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==1));...
    mean(cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==1));...
    mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==1) &...
    cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==1)) ;...
    mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==1))*...
    mean(cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==1));...

    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==1));...
    mean(cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==1));...

    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==1) &...
    cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==1)) ;...
    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==1))*...
    mean(cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==1))];




subplot(2,3,4);
bar([.9 1.1 1.4 1.6 2.4 2.6 2.9 3.1],bary');
set(gca,'XTick',[]), xlabel('CA1                    PFC');
ylabel(sprintf('Percentage of Task- \nResponsive Pyramidal Cells'))
% need to errorbar the grey bars (4 and 8)
% these can come from boostraps or just the bernoulli 95% likelihoods
% so binofit is the 95% conf, and binocdf is the likelihood
[~,pci]=binofit(round(bary(4)*sum(cellTable.isCA1==1 & cellTable.isPYR==1)),sum(cellTable.isCA1==1 & cellTable.isPYR==1));
[~,pci(2,:)]=binofit(round(bary(8)*sum(cellTable.isPFC==1 & cellTable.isPYR==1)),sum(cellTable.isPFC==1 & cellTable.isPYR==1));
pci=abs(pci-bary([4 8]));
hold on; errorbar([1.6 3.1],bary([4 8]),pci(:,1),pci(:,2),'k.','CapSize',0);
title('Beta and odorSelective overlap'); ylim([0 .4]);

% d. rr and selective overlap
bary=[mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==1));...
    mean(cellTable.respCoh(cellTable.isCA1==1 & cellTable.isPYR==1));...
    mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==1) &...
    cellTable.respCoh(cellTable.isCA1==1 & cellTable.isPYR==1)) ;...
    mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==1))*...
    mean(cellTable.respCoh(cellTable.isCA1==1 & cellTable.isPYR==1));...
   
    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==1));...
    mean(cellTable.respCoh(cellTable.isPFC==1 & cellTable.isPYR==1));... 
    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==1) &...
    cellTable.respCoh(cellTable.isPFC==1 & cellTable.isPYR==1)) ;...
    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==1))*...
    mean(cellTable.respCoh(cellTable.isPFC==1 & cellTable.isPYR==1))];


subplot(2,3,5);
bar([.9 1.1 1.4 1.6 2.4 2.6 2.9 3.1],bary');
set(gca,'XTick',[]), xlabel('CA1                    PFC');
ylabel(sprintf('Percentage of Task- \nResponsive Pyramidal Cells'))
% need to errorbar the grey bars (4 and 8)
% these can come from boostraps or just the bernoulli 95% likelihoods
% so binofit is the 95% conf, and binocdf is the likelihood
[~,pci]=binofit(round(bary(4)*sum(cellTable.isCA1==1 & cellTable.isPYR==1)),sum(cellTable.isCA1==1 & cellTable.isPYR==1));
[~,pci(2,:)]=binofit(round(bary(8)*sum(cellTable.isPFC==1 & cellTable.isPYR==1)),sum(cellTable.isPFC==1 & cellTable.isPYR==1));
pci=abs(pci-bary([4 8]));
hold on; errorbar([1.6 3.1],bary([4 8]),pci(:,1),pci(:,2),'k.','CapSize',0);
title('RR and odorSelective overlap');


% e. beta and choice INTERNEURONS

% d. rr and selective overlap
bary=[mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==0));...
    mean(cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==0));...
    mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==0) &...
    cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==0)) ;...
    mean(cellTable.odorP(cellTable.isCA1==1 & cellTable.isPYR==0))*...
    mean(cellTable.betaCoh(cellTable.isCA1==1 & cellTable.isPYR==0));...

    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==0));...
    mean(cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==0));...
    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==0) &...
    cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==0)) ;...
    mean(cellTable.odorP(cellTable.isPFC==1 & cellTable.isPYR==0))*...
    mean(cellTable.betaCoh(cellTable.isPFC==1 & cellTable.isPYR==0))];


subplot(2,3,6);
bar([.9 1.1 1.4 1.6 2.4 2.6 2.9 3.1],bary');
set(gca,'XTick',[]), xlabel('CA1                    PFC');
ylabel(sprintf('Percentage of Task- \nResponsive Interneurons'))
% need to errorbar the grey bars (4 and 8)
% these can come from boostraps or just the bernoulli 95% likelihoods
% so binofit is the 95% conf, and binocdf is the likelihood
[~,pci]=binofit(round(bary(4)*sum(cellTable.isCA1==1 & cellTable.isPYR==1)),sum(cellTable.isCA1==1 & cellTable.isPYR==1));
[~,pci(2,:)]=binofit(round(bary(8)*sum(cellTable.isPFC==1 & cellTable.isPYR==1)),sum(cellTable.isPFC==1 & cellTable.isPYR==1));
pci=abs(pci-bary([4 8]));
hold on; errorbar([1.6 3.1],bary([4 8]),pci(:,1),pci(:,2),'k.','CapSize',0);
title(sprintf('Beta and odorSelective overlap \n Interneurons'));


%% for pyrams vs ins
regions={'PFC','CA1'};
types={'pyr','in'};
respLockPct=table([0;0],[0;0],[0;0],[0;0],'VariableNames',{'PFCtot','PFClocked','CA1tot','CA1locked'},...
    'RowNames',types);
%contains({SuperRat(i).units.type},'in') & ...
for i=1:length(SuperRat)
    for j=1:length(regions)
        for k=1:length(types)
            mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) a(3)<.05, {SuperRat(i).units.taskResponsive}) &...
            contains({SuperRat(i).units.type},types{k}));
       
            respLockPct.([regions{j} 'tot'])(k)=respLockPct.([regions{j} 'tot'])(k)+length(mycells);
            respLockPct.([regions{j} 'locked'])(k)=respLockPct.([regions{j} 'locked'])(k)+...
                sum(cellfun(@(a) any(a<pcrit),{mycells.respRstat}));
        end
    end
end

openvar('respLockPct')
 %% grand average, CA1 all PFC all to local beta


regions={'PFC','CA1'};
types={'pyr','in'};
betaLockPct=table([0;0],[0;0],[0;0],[0;0],'VariableNames',{'PFCtot','PFClocked','CA1tot','CA1locked'},...
    'RowNames',types);
%contains({SuperRat(i).units.type},'in') & ...
for i=1:length(SuperRat)
    for j=1:length(regions)
        for k=1:length(types)
            mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) a(3)<.05, {SuperRat(i).units.taskResponsive}) &...
            contains({SuperRat(i).units.type},types{k}));
       
            betaLockPct.([regions{j} 'tot'])(k)=betaLockPct.([regions{j} 'tot'])(k)+length(mycells);
            respLobetaLockPctckPct.([regions{j} 'locked'])(k)=betaLockPct.([regions{j} 'locked'])(k)+...
                sum(cellfun(@(a) any(a<pcrit),{mycells.respRstat}));
        end
    end
end

openvar('betaLockPct')
%% panel E, number of beta and RR locked cells, pyr and IN on separate plots
% these are all to local (ca1 cells to ca1 lfp, pfc cells to pfc lfp)
regions={'PFC','CA1'}; 
cctype={'pyr','in'};
rhythm={'beta','resp'};



% now replicate claires analyses:

% first, she shows a few beta locked cells....
% This is a lenient threshold, so lets tabulate cells
pcrit=.05;
% top is pyr, bottom is IN
allcount=table([0;0],[0;0],[0;0],[0;0],[0;0],[0;0],'VariableNames',{'PFCtot','PFCbeta','PFCresp','CA1tot','CA1beta','CA1resp'},...
    'RowNames',{'pyr','in'});
%contains({SuperRat(i).units.type},'in') & ...

for i=1:length(SuperRat)
    for j=1:length(regions)
        % now for locating the table cells
        for k=1:2 % k is cell type
            mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) a(3)<pcrit, {SuperRat(i).units.taskResponsive}) & ...
            cellfun(@(a) contains(a,cctype{k}), {SuperRat(i).units.type}));       
        
            allcount.([regions{j} 'tot'])(k)=allcount.([regions{j} 'tot'])(k)+length(mycells);  % add this sess cells        
            for r=1:length(rhythm)
                allcount.([regions{j} rhythm{r}])(k)=allcount.([regions{j} rhythm{r}])(k)+sum(cellfun(@(a) a(j),{mycells.([rhythm{r} 'Rstat'])})<pcrit);
            end
        end
    end
end
 openvar('allcount')
 
 % now plot this as paired bars 
% so it will be ca1 pyr beta/resp, ca1 in beta/resp, pfc pyr beta/rsp
barpairs=[(allcount{:,[5 6]}./allcount{:,4})',(allcount{:,[2 3]}./allcount{:,1})' ];
 
hb=bar([0.9 1.1 1.4 1.6],barpairs',1,'LineStyle','none','FaceColor','flat');
hold on; plot([0.75 1.75],[pcrit pcrit],'k--')
for k=1:2, hb(k).CData=rhythmcolors(k,:); end
set(gca,'XTick',[0.9 1.1 1.4 1.6],'XTickLabel',{'pyr','in'});
xlabel('CA1                       PFC');
ylabel(sprintf('Proportion of \n Task Responsive Units')); box off;
legend('Beta','RR','chance','box', 'off');
xlim([.75 1.75]);
 
%% dial plots for phase preference of all cells across regions

% Panel C: polar plot of favored beta and rr phase in all task responsive neurons


% mean phase of all task responsive neuron
% task responsive is odor responsive (last row in odorresponsive)
rhythm={'beta','resp'};
regions={'PFC','CA1'};
for r=1:2
    figure;
cohtable=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCmean','PFCmvl','CA1mean','CA1mvl'});
pcrit=.05;
% this is all units and just for beta... we need to split pyr and in, and
% we need to run this for beta and RR

for i=1:length(SuperRat)
    for loc=1:length(regions) % for cell location
        % pull all cells in a location
        cellpool=SuperRat(i).units(contains({SuperRat(i).units.area},regions{loc}) &...
            cellfun(@(a) a(3)<pcrit,{SuperRat(i).units.taskResponsive}));
        % for LFP location
        if ~isempty(cellpool)
            for lfploc=1:length(regions)
                % get all task responsive cells
                respp=cellfun(@(a) a(3)<pcrit,{cellpool.taskResponsive});
                % get all cells with significant locking
                %betap=cellfun(@(a) a(betaloc),{cellpool.betaRstat});
                cohtable.([regions{loc} 'mean']){lfploc}=[cohtable.([regions{loc} 'mean']){lfploc}...
                    cellfun(@(a) a(lfploc),{cellpool(respp).([rhythm{r} 'mean'])})];
                cohtable.([regions{loc} 'mvl']){lfploc}=[cohtable.([regions{loc} 'mvl']){lfploc}...
                    cellfun(@(a) a(lfploc),{cellpool(respp).([rhythm{r} 'MVL'])})];
            end
        end
    end
end

% for each cell location and for each tetrode location

for i=1:2
    for k=1:2
        subplot(2,2,(i-1)*2+k);
        cohvals=cell2mat(cohtable.([regions{k} 'mean'])(i));
        p=polarhistogram(cohvals(~isnan(cohvals)),12,'Normalization','pdf');
        
        set(p,'LineStyle','none','FaceColor',colors(k,:));       
        hold on; set(gca,'ThetaTick',[0 90 180 270],'Rtick',[.1 .2 .3],'RTickLabel',{});
        [rho]=circ_mean(cohvals(~isnan(cohvals)));
        [mvl]=circ_r(cohvals(~isnan(cohvals)));
        [pval, z]=circ_rtest(cohvals(~isnan(cohvals)));
        PolarArrow(rho,mvl,[],colors(i,:).*.6);
        title(sprintf('Cells %s, %s %s \n n=%d mvl=%.2f p=%.2e',regions{k},rhythm{r}, regions{i},length(cohvals),z, pval));
    end
end
sgtitle(rhythm{r});
end
%% dial plots split out by pyrams and INs
% does beta coherence across cells group?

pcrit=.05;
useLocking=false;
for m=1:length(rhythm)
    for k=1:2
        cohtable=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCmean','PFCmvl','CA1mean','CA1mvl'},...
            'RowNames',{'PFCbeta','CA1beta'});
        for i=1:length(SuperRat)
            for loc=1:length(regions) % for cell location
                % pull all cells in a location
                cellpool=SuperRat(i).units(contains({SuperRat(i).units.area},regions{loc}) & ...
                    contains({SuperRat(i).units.type},types{k}));
                % for LFP location
                if ~isempty(cellpool)
                    for betaloc=1:length(regions)
                        % get all task responsive cells
                        okpool=cellfun(@(a) a(3)<pcrit,{cellpool.taskResponsive});

                        % get all cells with significant locking
                        if useLocking
                            okpool=okpool & cellfun(@(a) a(betaloc)<pcrit,{cellpool.([rhythm{m} 'Rstat'])});
                        end
                        cohtable.([regions{loc} 'mean']){betaloc}=[cohtable.([regions{loc} 'mean']){betaloc}...
                            cellfun(@(a) a(betaloc),{cellpool(okpool).([rhythm{m} 'mean'])})];
                        cohtable.([regions{loc} 'mvl']){betaloc}=[cohtable.([regions{loc} 'mvl']){betaloc}...
                            cellfun(@(a) a(betaloc),{cellpool(okpool).([rhythm{m} 'MVL'])})];
                    end
                end
            end
        end

        % for each cell location and for each tetrode location
        figure;
        for i=1:2
            for j=1:2
                subplot(2,2,(i-1)*2+j);
                betavals=cell2mat(cohtable.([regions{j} 'mean'])(i));
                p=polarhistogram(betavals(~isnan(betavals)),12,'Normalization','pdf');
                set(p,'LineStyle','none','FaceColor',colors(j,:));
                hold on;
                [rho]=circ_mean(betavals(~isnan(betavals)));
                [mvl]=circ_r(betavals(~isnan(betavals)));
                PolarArrow(rho,mvl,[],colors(i,:).*.6);
                pval=circ_rtest(betavals(~isnan(betavals)));
                set(gca,'ThetaTick',[0 90 180 270],'Rtick',[.1 .2 .3],'RTickLabel',{});

                title(sprintf('Cells %s, %s %s \n n=%d, p=%.2e',regions{j}, rhythm{m}, regions{i},length(betavals),pval));

            end
        end
        sgtitle(sprintf('%s Coherence %s',rhythm{m}, types{k}));
    end
end



%%
% Panel C: polar plot of favored beta and rr phase in all task responsive
% neurons AND COHERENT CELLS

% PYR AND IN BREAKOUT

% mean phase of all task responsive neuron
% task responsive is odor responsive (last row in odorresponsive)
rhythm={'beta','resp'};
regions={'PFC','CA1'};
pcrit=.05;
for ct=1:2
    for r=1:2
        figure;
        cohtable=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCmean','PFCmvl','CA1mean','CA1mvl'});

        % this is all units and just for beta... we need to split pyr and in, and
        % we need to run this for beta and RR
        
        for i=1:length(SuperRat)
            for loc=1:length(regions) % for cell location
                % pull all cells in a location
                cellpool=SuperRat(i).units(contains({SuperRat(i).units.area},regions{loc}) &...
                    contains({SuperRat(i).units.type},types{ct}));
                % for LFP location
                if ~isempty(cellpool)
                    for lfploc=1:length(regions)
                        % get all task responsive cells
                        respp=cellfun(@(a) a(3)<.05,{cellpool.taskResponsive});% &...
                            %cell2mat(cellfun(@(a) a(lfploc)<pcrit, {cellpool.([rhythm{r} 'Rstat'])},'UniformOutput',false));
                        % get all cells with significant locking
                        %betap=cellfun(@(a) a(betaloc),{cellpool.betaRstat});
                        cohtable.([regions{loc} 'mean']){lfploc}=[cohtable.([regions{loc} 'mean']){lfploc}...
                            cellfun(@(a) a(lfploc),{cellpool(respp).([rhythm{r} 'mean'])})];
                        cohtable.([regions{loc} 'mvl']){lfploc}=[cohtable.([regions{loc} 'mvl']){lfploc}...
                            cellfun(@(a) a(lfploc),{cellpool(respp).([rhythm{r} 'MVL'])})];
                    end
                end
            end
        end
        
        % for each cell location and for each tetrode location
        
        for i=1:2
            for k=1:2
                subplot(2,2,(i-1)*2+k);
                cohvals=cell2mat(cohtable.([regions{k} 'mean'])(i));
                p=polarhistogram(cohvals(~isnan(cohvals)),12,'Normalization','pdf');
                pval=circ_rtest(cohvals(~isnan(cohvals)));
                title(sprintf('Cells %s, %s %s \n n=%d, p=%.2e',...
                    regions{k},rhythm{r}, regions{i},sum(~isnan(cohvals)), pval));
                set(p,'LineStyle','none','FaceColor',colors(k,:));
                hold on; set(gca,'ThetaTick',[0 90 180 270],'Rtick',[.1 .2 .3],'RTickLabel',{});
                [rho]=circ_mean(cohvals(~isnan(cohvals)));
                [mvl]=circ_r(cohvals(~isnan(cohvals)));
                PolarArrow(rho,mvl,[],colors(i,:).*.6);
                
            end
        end
        sgtitle(sprintf('%s %s',rhythm{r},types{ct}));
    end
end
%%
%
%
%
%
%
%     correct vs incorrect
%
%
%
%


%% Beta locking higher for correct vs incorrect?

% and option for only locking


%figure; subplot(2,1,1); histogram(bootC); subplot(2,1,2); histogram(bootI);
%linkaxes(get(gcf,'Children'),'x');
%{
 now tabulate all the cells first get hpc to hpc, pfc pfc

changing the pcrit for the task responsiveness or the coherence metric does
not meaningfully change the main result.



%}
onlysig=1; saveout=0;
pcrit=.05;
if saveout, savefolder=uigetdir; end

celltype={'pyr','in'}; 
regions={'PFC','CA1'};
LFPregs={'PFC','CA1','OB'};
rhythm={'beta','resp'}; 

%%
clear allyy pval
for r=1:2 % rhythm
    for ct=1:2  % celltype
        allcount=table({[];[];[]},{[];[];[]},{[];[];[]},{[];[];[]},'VariableNames',{'PFCvals','PFCdprimes','CA1vals','CA1dprimes'});
        for i=1:length(SuperRat) % session
            for j=1:length(regions) % region for units
                for k=1:3 % region for LFP
                    % just take odor responsive cells in that regin that cohere to beta
                    if onlysig==1
                        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
                            cellfun(@(a) a(k)<=pcrit, {SuperRat(i).units.([lower(rhythm{r}) 'Rstat'])}) &...
                            cellfun(@(a) a(3)<=.05, {SuperRat(i).units.taskResponsive}) &...
                            contains({SuperRat(i).units.type},celltype{ct}));
                    else
                        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
                            cellfun(@(a) a(3)<=.05, {SuperRat(i).units.taskResponsive}) &...
                            contains({SuperRat(i).units.type},celltype{ct}));
                    end
                    if ~isempty(mycells)
                        sessvals=cell2mat(cellfun(@(a) a(k,:), {mycells.([rhythm{r} 'PhaseDprime'])},'UniformOutput',0)');
                        allcount.([regions{j} 'vals']){k}=[allcount.([regions{j} 'vals']){k} ; sessvals(:,1:2)];
                        allcount.([regions{j} 'dprimes']){k}=[allcount.([regions{j} 'dprimes']){k}; sessvals(:,1)-sessvals(:,2)];
                    end
                end
            end
        end
        
        for i=1:2 % cell region (cols)
            figure('position',[1100,450-80*((i-1)*2+k),600 450]); sp=[]; sp2=[];
            for k=1:3 % LFP region (rows)
                try

                    sp(k)=subplot(2,3,k);
                    %bar([1 2],mean(allcount.([regions{i} 'vals']){k},'omitnan'),'LineStyle','none','FaceColor',colors(i,:)); hold on;
                    plot(repmat([1;2],1,length(allcount.([regions{i} 'vals']){k})),allcount.([regions{i} 'vals']){k}',...
                        '-','color',rhythmcolors(r,:),'LineWidth',3);
                    ylabel(sprintf('Rate Adjusted \n Vector Length'));
                    box off;

                    xlim([0.5 2.5]); set(gca,'XTick',[1 2],'XTickLabel',{'Correct','Incorrect'});
                    sp2(k)=subplot(2,3,k+3); [y,x]=histcounts(allcount.([regions{i} 'dprimes']){k},-.3:.04:.3,'Normalization','probability');
                    barh(x(2:end),y,'LineStyle','none','BarWidth',1,'FaceColor',colors(i,:));
                    box off; hold on; plot(get(gca,'XLim'),[0 0],'k'); %repmat(mean(allcount.([regions{i} 'dprimes']){k},'omitnan'),1,2),'k');
                    xlabel('Proportion of units'); ylabel('Change in MVL'); box off;
                    [p]=signrank(diff(allcount.([regions{i} 'vals']){k}'));
                    mymeans=mean(allcount.([regions{i} 'vals']){k},'omitnan');
                    title(sprintf('%s type:%s n=%d,\n %s LFP in %s \n C=%.2f I=%.2f, \n signrank=%.2e',...
                        regions{i},celltype{ct},length((allcount.([regions{i} 'vals']){k})),...
                        rhythm{r}, LFPregs{k},mymeans(1),mymeans(2),p));

                    fprintf('%s type:%s n=%d, %s LFP in %s \n C=%.2f I=%.2f,signrank=%.2e\n',...
                        regions{i},celltype{ct},length((allcount.([regions{i} 'vals']){k})),...
                        rhythm{r}, LFPregs{k},mymeans(1),mymeans(2),p);

                end
                
            end
            ylim([-.35 .35]);
            subplot(2,3,1); ylim([0 .8]);
            linkaxes(sp); linkaxes(sp2);
            fprintf('\n'); drawnow;
            if saveout
                % make sure you track whether these are coherent
                % cells or not
                if onlysig
                    figname=fullfile(savefolder,sprintf('CorrectIncorrect_%sCohCells_%s_from_%s',rhythm{r},celltype{ct},regions{i}));
                else
                    figname=fullfile(savefolder,sprintf('CorrectIncorrect_taskRespCells%s_%s_from_%s',rhythm{r},celltype{ct},regions{i}));
                end
                print(figname,'-dsvg');
                kill;
            end
            % end figure/region
        end
        % end celltype
    end
    % end rhythm
end

%% which regions beta do cells cohere best to?
allCells=cell2mat({SuperRat.units});
okcells=cell2mat(cellfun(@(a) a, {allCells.taskResponsive}, 'UniformOutput', false)');
allCells=allCells(okcells(:,3)<.05); % only take task responsive cells
% lets start by plot3 scatter
%
pcrit=.05/3;
figure;
for ct=1:2 % celltype
    for cellReg=1:2 % source region
        mycells=allCells(contains({allCells.area},regions{cellReg}) &...
            cellfun(@(a) any(a<=pcrit), {allCells.betaRstat}) &...
            contains({allCells.type},types{ct}));
        % what index is the max mvl?
        allMVL=cell2mat({mycells.betaRstat}');
        [~,favind]=min(allMVL,[],2);
        subplot(2,2,(ct-1)*2+cellReg);
        pie(histcounts(favind,[.5:1:3.5])); 
        title(sprintf('%s units in %s \n n=%d',types{ct},regions{cellReg},length(favind)))
    end
end
legend({'PFCbeta','CA1beta','OBbeta'});
% and for rr
figure;
for ct=1:2 % celltype
    for cellReg=1:2 % source region
        mycells=allCells(contains({allCells.area},regions{cellReg}) &...
            cellfun(@(a) any(a<=pcrit), {allCells.respRstat}) &...
            contains({allCells.type},types{ct}));
        % what index is the max mvl?
        allMVL=cell2mat({mycells.respRstat}');
        [~,favind]=min(allMVL,[],2);
        subplot(2,2,(ct-1)*2+cellReg);
        pie(histcounts(favind,[.5:1:3.5])); 
        title(sprintf('%s units in %s \n n=%d',types{ct},regions{cellReg},length(favind)))
    end
end
legend({'PFCrr','CA1rr','OBrr'});
%INcohN=INcoh-mean(INcoh,2);

%% and now pairwise
% it'll be a four row by three column
orders=[1 2; 1 3; 2 3];
figure
for ct=1:2 % celltype
    for cellReg=1:2 % source region
        mycells=allcells(contains({allCells.area},regions{cellReg}) &...
            cellfun(@(a) any(a<=pcrit), {allCells.betaRstat}) &...
            contains({allCells.type},types{ct}));
        allMVL=cell2mat({mycells.betaMVL}');
        means=mean(allMVL); errs=SEM(allMVL,1);
        for pr=1:length(orders)
            subplot(4,3,((ct-1)*2+cellReg-1)*3+pr);
            scatter(allMVL(:,orders(pr,1)),allMVL(:,orders(pr,2)),12,'filled');
            lims=max([get(gca,'YLim') get(gca,'XLim')]);
            hold on; plot([0 max(lims)],[0 max(lims)],'k--');
            % and now the mean value target
            plot(means([orders(pr,1) orders(pr,1)]),(means(orders(pr,2))+[errs(orders(pr,2)), -errs(orders(pr,2))]),...
                'color','r')
            plot((means(orders(pr,1))+[errs(orders(pr,1)), -errs(orders(pr,1))]),...
                means([orders(pr,2) orders(pr,2)]),'color','r')
            title(sprintf('%s units in %s \n p=%.2e',types{ct},regions{cellReg},...
                signrank(allMVL(:,orders(pr,1))-allMVL(:,orders(pr,2)))));
            xlabel(regions2{orders(pr,1)}); ylabel(regions2{orders(pr,2)});
            % grab coding depth for first pair
        end
    end
end
sgtitle('beta')

figure;
for ct=1:2 % celltype
    for cellReg=1:2 % source region
        mycells=allcells(contains({allCells.area},regions{cellReg}) &...
            cellfun(@(a) any(a<=pcrit), {allCells.respRstat}) &...
            contains({allCells.type},types{ct}));
        allMVL=cell2mat({mycells.respMVL}');
        means=mean(allMVL); errs=SEM(allMVL,1);
        for pr=1:length(orders)
            subplot(4,3,((ct-1)*2+cellReg-1)*3+pr);
            scatter(allMVL(:,orders(pr,1)),allMVL(:,orders(pr,2)),12,'filled');
            lims=max([get(gca,'YLim') get(gca,'XLim')]);
            hold on; plot([0 max(lims)],[0 max(lims)],'k--');
            % and now the mean value target
            plot(means([orders(pr,1) orders(pr,1)]),(means(orders(pr,2))+[errs(orders(pr,2)), -errs(orders(pr,2))]),...
                'color','r')
            plot((means(orders(pr,1))+[errs(orders(pr,1)), -errs(orders(pr,1))]),...
                means([orders(pr,2) orders(pr,2)]),'color','r')
            title(sprintf('%s units in %s \n p=%.2e',types{ct},regions{cellReg},...
                signrank(allMVL(:,orders(pr,1))-allMVL(:,orders(pr,2)))));
            xlabel(regions2{orders(pr,1)}); ylabel(regions2{orders(pr,2)});
            % grab coding depth for first pair
        end
    end
end
sgtitle('respiratory rhythm')

%%
%
%
%
%
%
%      lara rangel coherence metrics, unfinished
%
%
%
%
%
%% the lara rangel way, using coherencycpt from chronux


addpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\MyChronux'));
for i=1:length(SuperRat)
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    corrIncorr=SuperRat(i).trialdata.CorrIncorr10;

    for j=1:length(SuperRat(i).units)
        % get this cells spikes
        spikeTimes=SuperRat(i).units(j).ts(:,1);
        % cell mat of spiketimes
        for k=1:length(regions)

            % get Beta from region K
            Betadata=SuperRat(i).([regions{k} 'beta']);
            % can cellfun this- pull the filtered data
            params=struct('tapers',[3 5],'fs',1500,'pad',1);
            [C,phi,S1,S2,S12] = calcSpikeFieldCoherence(spikeTimes,Betadata(:,[1 2]),trialTimes,params);
            % this is for all the trials,
            
            
            % here we can do it with downsampling
            
            % now for each task responsive cell, grab trialPhases.
            % then we can bootstrap a mean phase, calculate the difference
            % in phases, and then bootstrap a difference in MVL
            
            % save out for each cell a 
 
        end
    end
end


rmpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\MyChronux'));
%%
% calculating when beta becomes good on each trial

% i think the gist of this is that the beta power just keeps going up the
% longer the animal is in the odor port...


for i=1:length(SuperRat)
    % get odor starts
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    odorTimes([odorTimes(:,2)-odorTimes(:,1)]<.5,:)=[];
    figure;
    for j=1:length(regions)
        betaPow=SuperRat(i).([regions{j} 'beta']);
        betaPow(:,4)=zscore(betaPow(:,4));
        respPow=SuperRat(i).([regions{j} 'resp']);
        respPow(:,4)=zscore(respPow(:,4));
        
        % now gather the amplitudes
        [~,~,~,lfpInds,~,times]=event_spikes(betaPow(:,1),odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
        betaData=cellfun(@(a) betaPow(a,4), lfpInds,'UniformOutput',false);
        respData=cellfun(@(a) respPow(a,4), lfpInds,'UniformOutput',false);
        
        subplot(4,3,j)
        for k=1:length(betaData)
            plot(times{k},betaData{k},'k');
            hold on;
        end
        title(sprintf('%s Beta Power locked to sniff start',regions{j}));
        subplot(4,3,j+3)
        for k=1:length(respData)
            plot(times{k},respData{k},'k');
            hold on;
        end
        title(sprintf('%s RR Power locked to sniff start',regions{j}));
        sgtitle(sprintf('%d rat %s ses %d',i,SuperRat(i).name,SuperRat(i).daynum));
        
        subplot(4,3,j+6)
        durs=cell2mat(cellfun(@(a) length(a), betaData,'UniformOutput',false))/1000;
        maxxes=cell2mat(cellfun(@(a) max([nan find(a>2,1,'first')]), betaData,'UniformOutput',false))/1000;
        scatter(durs,maxxes); title('beta max vs sniffdur');
        
        subplot(4,3,j+9)
        durs=cell2mat(cellfun(@(a) length(a), respData,'UniformOutput',false))/1000;
        maxxes=cell2mat(cellfun(@(a) max([nan find(a>2,1,'first')]), respData,'UniformOutput',false))/1000;
        scatter(durs,maxxes); title('resp max vs sniffdur');
        linkaxes(get(gcf,'Children'));
    end
end

    
%% so we know that beta power increases with time, but what about spike-beta coherence
% all ca1 cells tend to fire at the same beta phase, so maybe the better
% clustering leads to quicker decisions

% this isnt tractable to test though because you need an instantaneous
% estimate of phase clustering


% lets start with pyramdial cells
celltype='pyr';

for i=1:length(SuperRat)
    % first lets focus on local beta.  To do this we'll do each region
    % separately
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    odorTimes([odorTimes(:,2)-odorTimes(:,1)]<.5,:)=[];
    %figure;
    for j=1:length(regions)
        betaData=SuperRat(i).([regions{j} 'beta']);
        betaData(:,4)=zscore(betaData(:,4));
        respData=SuperRat(i).([regions{j} 'resp']);
        respData(:,4)=zscore(respData(:,4));
        
        mycells=SuperRat(i).units(cellfun(@(a) contains(a,regions{j}), {SuperRat(i).units.area}) &...
            cellfun(@(a) contains(a,celltype), {SuperRat(i).units.type}) & ...
            cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
        
        allSpikes=sort(cell2mat({mycells.ts}'));
        [~,~,~,~,odorSpikes,evSpikes]=event_spikes(allSpikes,odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
        
        % now grab the phases
        spikePhases=cellfun(@(a) interp1(betaData(:,1),betaData(:,3),a,'nearest'), odorSpikes,'UniformOutput',false);
        
        % now find when the mvl reaches significance... maybe histogram
        % them?
        % for each trial, divide into 20 spikes, get mean time of those
        % spikes, and calc r.  find first that r is sig
        phaseTimes=nan(1,length(spikePhases));
        for q=1:length(spikePhases)
            bins=1:10:length(spikePhases{q});
            binmean=[]; binr=[];
            for r=1:length(bins)-1
                    binmean(r)=mean(evSpikes{q}(bins(r):bins(r+1)));
                    binr(r)=circ_rtest(spikePhases{q}(bins(r):bins(r+1)));
                    if binr<.05, break; end
            end
            firstind=find(binr<.05,1,'first');
            if ~isempty(firstind), phaseTimes(q)=binmean(firstind);
            else; phaseTimes(q)=nan;
            end
        end
        
        figure; scatter(odorTimes(:,2)-odorTimes(:,1),phaseTimes);
    end
end
        

%% lets confirm we have aligned data.  First lets look at Beta-odor period alignment
% this just plots the beta amplitude locked to sniff start timestamps for a
% few sessions

%{
for i=5:8 %length(SuperRat)
    % get sniffstarts
    figure;
    sniffend=SuperRat(i).trialdata.sniffend;
    lfpAmp=zscore(SuperRat(i).PFCbeta(:,4));
    betaim=[];
    for k=1:length(sniffend)
        betastart=find(SuperRat(i).PFCbeta(:,1)>sniffend(k),1,'first');
        betaim(k,:)=lfpAmp(betastart-3000:betastart+3000);
    end
    subplot(1,2,1);
    imagesc(-2:(1/1500):2,length(sniffend),betaim);
    yyaxis right; plot(-2:(1/1500):2,mean(betaim),'k','LineWidth',2);
    title(sprintf('%s day %d PFC',SuperRat(i).name,SuperRat(i).daynum));
    lfpAmp=zscore(SuperRat(i).CA1beta(:,4));
    betaim=[];
    for k=1:length(sniffend)
        betastart=find(SuperRat(i).PFCbeta(:,1)>sniffend(k),1,'first');
        betaim(k,:)=lfpAmp(betastart-3000:betastart+3000);
    end
    subplot(1,2,2);
    imagesc(-2:(1/1500):2,length(sniffend),betaim);
    yyaxis right; plot(-2:(1/1500):2,mean(betaim),'k','LineWidth',2);
    title(sprintf('%s day %d CA1',SuperRat(i).name,SuperRat(i).daynum));
end
    
%}

    