% ClaireRRanalysis



has not been edited to use rr instead of beta yet!!!!!!!!!!!



%% Comments block
%{


Claires algorithm:
1. find tet with most cells
2. get that tets RR filtered lfp
3. interp eeg times

4. get spikes in odor delivery period ONLY
        or epoch spikes based on when beta envelope is >2sd above mean,
        using a fwhm approach.
5. get phase of those spikes
6. calc mean phase, mean vector length, p value
7. bootstrap that dataset (could probably pull beta from different trials)

8. aggregate for pyrams and INs and pfc and hpc





%}
regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')];

%end
%% now run the stats on the cells

%

for i=1:length(SuperRat)
    myclock=tic;
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorr=SuperRat(i).trialdata.CorrIncorr10;
    
    for j=1:length(SuperRat(i).units)
        myspikes=EpochCoords(SuperRat(i).units(j).ts(:,1),snifftimes(trialCorr==1,:));

        for k=1:length(regions)
            myLFP=SuperRat(i).([regions{k} 'resp']);
            % insert here case if no viable LFP
            SuperRat(i).units(j).respmean(k)=nan;
            SuperRat(i).units(j).respMVL(k)=nan;
            SuperRat(i).units(j).respOstat(k)=nan;
            SuperRat(i).units(j).respRstat(k)=nan;
            
            if length(myspikes)>5 && ~isempty(myLFP)
                
                % if there are enough spikes to analyze, get the mean
                % phase, the mvl, the r test and O test
                phasesCorr=interp1(myLFP(:,1),myLFP(:,3),myspikes,'nearest');
                SuperRat(i).units(j).respmean(k)=circ_mean(phasesCorr);
                SuperRat(i).units(j).respMVL(k)=circ_r(phasesCorr);
                SuperRat(i).units(j).respOstat(k)=circ_otest(phasesCorr);
                SuperRat(i).units(j).respRstat(k)=circ_rtest(phasesCorr);

                % plot it out???
%                 
%                 if SuperRat(i).units(j).respRstat(k)<.01
%                     figure; ha=histogram([phasesCorr; phasesCorr+pi*2],linspace(-pi,pi*3,10+ceil(sqrt(length(phasesCorr)))),...
%                         'Normalization','probability');
%                     set(ha,'LineStyle','none','FaceColor',colors(k,:));
%                     set(gca,'XTick',[pi.*[-1:4]],'XTickLabel',{'-2\pi','','0','','2\pi'});
%                     ylabel('Probability of Spike'); xlabel('Beta Phase'); box off;
%                      title(sprintf(' %s cell in %s, lfp %s \n n=%.f, mean %.2f, mvl=%.2f, p=%.2e',...
%                         SuperRat(i).units(j).type, SuperRat(i).units(j).area, regions{k}, length(phasesCorr),...
%                         SuperRat(i).units(j).respmean(k),SuperRat(i).units(j).respMVL(k),SuperRat(i).units(j).respRstat(k)));
%                     figure;
%                     scatter(repmat(c,2,1),[phasesCorr; phasesCorr+pi*2]);
%                     box off;
%                 end
            end
        end
    end
end


%%

%
%
%
%  Correct vs incorrect
%
%
%



%% now measuring the strength of this modulation.

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
nboots=200; bootclock=tic;
for i=1:length(SuperRat)
    myclock=tic;
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorr=SuperRat(i).trialdata.CorrIncorr10;
    
    h=waitbar(0,sprintf('Starting animal %s %d', SuperRat(i).name,SuperRat(i).daynum));
    for j=1:length(SuperRat(i).units)
        [~,~,~,~,myspikes]=event_spikes(SuperRat(i).units(j).ts(:,1),snifftimes(:,1),0,snifftimes(:,2)-snifftimes(:,1));
        
        for k=1:length(regions)
            myLFP=SuperRat(i).([regions{k} 'resp']);
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


%%
%
%
%
%
%        Plot functions
%
%
%
%
%
% what we want to plot is...
% 1. the % change in coherence across cells from correct to incorrect-
% something like signed difference over sum (to adjust for overall
% coherence)
% 2. the %% of cells whose coherence significantly chagnes in one direction


% I think the result is likely that ca1 INS change their coherence the
% most, and that this is caused by PFC top-down connections that then bias
% the behavior of the animal.  So the question is can we even see the
% mechanism by which the bottom up interaction occurs?
%% how many coherent cells in each region?
% now replicate claires analyses:
% local coherence and remote (cell in region 1, lfp region2)

% first, she shows a few beta locked cells....
% This is a lenient threshold, so lets tabulate cells
regions={'PFC','CA1'};
allcount=table([0;0],[0;0],[0;0],[0;0],'VariableNames',{'PFCtot','PFClocked','CA1tot','CA1locked'});
%contains({SuperRat(i).units.type},'in') & ...

for i=1:length(SuperRat)
    for j=1:length(regions)
        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
        for k=1:length(regions)
            allcount.([regions{j} 'tot'])(k)=allcount.([regions{j} 'tot'])(k)+length(mycells);
            allcount.([regions{j} 'locked'])(k)=allcount.([regions{j} 'locked'])(k)+sum(cellfun(@(a) a(k),{mycells.respRstat})<.05);
        end
    end
end
openvar('allcount')
%% for pyrams vs ins
regions={'PFC','CA1'};
types={'pyr','in'};
allcount=table([0;0],[0;0],[0;0],[0;0],'VariableNames',{'PFCtot','PFClocked','CA1tot','CA1locked'});
%contains({SuperRat(i).units.type},'in') & ...
for i=1:length(SuperRat)
    for j=1:length(regions)
        for k=1:length(types)
        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}) &...
            contains({SuperRat(i).units.type},types{k}));
       
            allcount.([regions{j} 'tot'])(k)=allcount.([regions{j} 'tot'])(k)+length(mycells);
            allcount.([regions{j} 'locked'])(k)=allcount.([regions{j} 'locked'])(k)+sum(cellfun(@(a) a(k),{mycells.respRstat})<.05);
        end
    end
end
openvar('allcount')

%% 
% Panel C:
% mean phase of all task responsive neuron
% task responsive is odor responsive (last row in odorresponsive)

betatable=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCmean','PFCmvl','CA1mean','CA1mvl'},...
    'RowNames',{'PFCbeta','CA1beta'});
pcrit=.05;
for k=1:2

for i=1:length(SuperRat)
    for loc=1:length(regions) % for cell location
        % pull all cells in a location
        cellpool=SuperRat(i).units(contains({SuperRat(i).units.area},regions{loc}) & ...
            contains({SuperRat(i).units.type},types{k}));
        % for LFP location
        if ~isempty(cellpool)
            for betaloc=1:length(regions)
                % get all task responsive cells
                betap=cellfun(@(a) any(a(:,4)<pcrit),{cellpool.OdorResponsive});
                
                
                % get all cells with significant locking
                %betap=cellfun(@(a) a(betaloc),{cellpool.betaRstat});
                betatable.([regions{loc} 'mean']){betaloc}=[betatable.([regions{loc} 'mean']){betaloc}...
                    cellfun(@(a) a(betaloc),{cellpool(betap).betamean})];
                betatable.([regions{loc} 'mvl']){betaloc}=[betatable.([regions{loc} 'mvl']){betaloc}...
                    cellfun(@(a) a(betaloc),{cellpool(betap).betaMVL})];
            end
        end
    end
end

% for each cell location and for each tetrode location
figure;
for i=1:2
    for j=1:2
        subplot(2,2,(i-1)*2+j);
        betavals=cell2mat(betatable.([regions{j} 'mean'])(i));
        p=polarhistogram(betavals(~isnan(betavals)),12,'Normalization','pdf');
        [pval]=circ_rtest(betavals);
        title(sprintf('Cells %s, beta %s \n p=%.2e',regions{j},regions{i}));
        set(p,'LineStyle','none','FaceColor',colors(j,:));       
        hold on;
        [rho]=circ_mean(betavals(~isnan(betavals)));
        [mvl]=circ_r(betavals(~isnan(betavals)));
        PolarArrow(rho,mvl,[],colors(i,:).*.6);
        
    end
end
sgtitle(sprintf('%s',types{k}));
end
%% my panel c looks way different from hers

% panel E.

% can only do Beta here


allcount=table(0,0,0,0,'VariableNames',{'PFCtot','PFClocked','CA1tot','CA1locked'});
pcrit=.05;

for i=1:length(SuperRat)
    for j=1:length(regions)
        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}));
        
        allcount.([regions{j} 'tot'])=allcount.([regions{j} 'tot'])+length(mycells);
        coders=cellfun(@(a) any(a<=pcrit),{mycells.betaRstat});
        allcount.([regions{j} 'locked'])=allcount.([regions{j} 'locked'])+sum(coders);
    end
end
figure;
bar([ allcount.CA1locked/allcount.CA1tot allcount.PFClocked/allcount.PFCtot]);
set(gca,'XTickLabel',{'PFC','CA1'});
ylabel('Percent Beta Coherent'); xlabel('Cell & Beta source');
ylim([0 1]);
%% This is for all task responsive cells

%figure; subplot(2,1,1); histogram(bootC); subplot(2,1,2); histogram(bootI);
%linkaxes(get(gcf,'Children'),'x');

% now tabulate all the cells
% first get hpc to hpc, pfc pfc

%savefolder=uigetdir();
allcount=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCvals','PFCdprimes','CA1vals','CA1dprimes'});
celltype='in';

for i=1:length(SuperRat)
    for j=1:length(regions)
        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) any(a(:,4)<.05), {SuperRat(i).units.OdorResponsive}) &...
            contains({SuperRat(i).units.type},celltype));
        if ~isempty(mycells)
        for k=1:length(regions)
            sessvals=cell2mat(cellfun(@(a) a(k,:), {mycells.betaPhaseDprime},'UniformOutput',0)');
            allcount.([regions{j} 'vals']){k}=[allcount.([regions{j} 'vals']){k} ; sessvals(:,1:2)];
            allcount.([regions{j} 'dprimes']){k}=[allcount.([regions{j} 'dprimes']){k}; sessvals(:,1)-sessvals(:,2)];
        end
        end
    end
end

for i=1:2
    for k=1:2

        figure('position',[1100,1050, 500 300]);  subplot(1,2,1);
        bar([1 2],mean(allcount.([regions{i} 'vals']){k},'omitnan'),'LineStyle','none','FaceColor',colors(i,:)); hold on;
        plot(repmat([1;2],1,length(allcount.([regions{i} 'vals']){k})),allcount.([regions{i} 'vals']){k}','o-','color',[.7 .7 .7]);
        ylabel(sprintf('Rate Adjusted \n Vector Length'));
        box off
        xlim([0.5 2.5]); set(gca,'XTickLabel',{'Correct','Incorrect'});
        subplot(1,2,2); [y,x]=histcounts(allcount.([regions{i} 'dprimes']){k},'Normalization','probability');
        barh(x(2:end),y,'LineStyle','none','BarWidth',1,'FaceColor',colors(i,:));
        hold on; plot(get(gca,'XLim'),[0 0],'k'); %repmat(mean(allcount.([regions{i} 'dprimes']){k},'omitnan'),1,2),'k');
        xlabel('Proportion of units'); ylabel('Change in MVL'); box off;
        [p]=signrank(diff(allcount.([regions{i} 'vals']){k}'));
        mymeans=mean(allcount.([regions{i} 'vals']){k},'omitnan');
        sgtitle(sprintf('All Task %s cells in %s LFP in %s \n C=%.2f I=%.2f,signrank=%.2e',...
            celltype,regions{i},regions{k},mymeans(1),mymeans(2),p));
        %figname=fullfile(savefolder,sprintf('BetaCoh All Task %s cells from %s LFP in %s',celltype,regions{i},regions{k}));
        %print(figname,'-dsvg')
    end
    
end
% figure will be something like
% first paired points, C, IC, then a histogram beside it for difference
% so I need concatenated HPC-HPC C, IC, Dprime


%% This is for just Beta Locking cells (pval<.05)


%figure; subplot(2,1,1); histogram(bootC); subplot(2,1,2); histogram(bootI);
%linkaxes(get(gcf,'Children'),'x');
%{
 now tabulate all the cells first get hpc to hpc, pfc pfc

changing the pcrit for the task responsiveness or the coherence metric does
not meaningfully change the main result.



%}
pcrit=.05;
celltype={'pyr','in'};
for ct=1:2
    allcount=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCvals','PFCdprimes','CA1vals','CA1dprimes'});
    for i=1:length(SuperRat)
        for j=1:length(regions)
            for k=1:length(regions)
                % just take odor responsive cells in that regin that cohere to beta
                mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
                    cellfun(@(a) any(a(:,4)<=pcrit), {SuperRat(i).units.OdorResponsive}) &...
                    cellfun(@(a) a(k)<=pcrit, {SuperRat(i).units.betaRstat}) &...
                    contains({SuperRat(i).units.type},celltype{ct}));
                if ~isempty(mycells)
                    sessvals=cell2mat(cellfun(@(a) a(k,:), {mycells.BetaPhaseDprime},'UniformOutput',0)');
                    allcount.([regions{j} 'vals']){k}=[allcount.([regions{j} 'vals']){k} ; sessvals(:,1:2)];
                    allcount.([regions{j} 'dprimes']){k}=[allcount.([regions{j} 'dprimes']){k}; sessvals(:,1)-sessvals(:,2)];
                end
            end
        end
    end
    
    for i=1:2
        for k=1:2
            figure('position',[1100,1050-80*((i-1)*2+k), 500 300]);
            subplot(1,2,1);
            bar([1 2],mean(allcount.([regions{i} 'vals']){k},'omitnan'),'LineStyle','none','FaceColor',colors(i,:)); hold on;
            plot(repmat([1;2],1,length(allcount.([regions{i} 'vals']){k})),allcount.([regions{i} 'vals']){k}','o-','color',[.7 .7 .7]);
            ylabel(sprintf('Rate Adjusted \n Vector Length'));
            box off;
            
            xlim([0.5 2.5]); set(gca,'XTickLabel',{'Correct','Incorrect'});
            subplot(1,2,2); [y,x]=histcounts(allcount.([regions{i} 'dprimes']){k},'Normalization','probability');
            barh(x(2:end),y,'LineStyle','none','BarWidth',1,'FaceColor',colors(i,:));
            box off; hold on; plot(get(gca,'XLim'),[0 0],'k'); %repmat(mean(allcount.([regions{i} 'dprimes']){k},'omitnan'),1,2),'k');
            xlabel('Proportion of units'); ylabel('Change in MVL'); box off;
            [p]=signrank(diff(allcount.([regions{i} 'vals']){k}'));
            mymeans=mean(allcount.([regions{i} 'vals']){k},'omitnan');
            title(sprintf('B coh Task %s cells in %s LFP in %s \n C=%.2f I=%.2f,signrank=%.2e',...
                celltype{ct},regions{i},regions{k},mymeans(1),mymeans(2),p));
            %figname=fullfile(savefolder,sprintf('BetaCoh Task and Coherent %s cells from %s LFP in %s',celltype,regions{i},regions{k}));
            %print(figname,'-dsvg')
            
        end
    end
end
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
    