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

% set region parameters here!!
regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')];
types={'pyr','in'};

%% Overall beta coherence using whole odor period and just the correct trials
% we should be using stats only on local beta, but we can do all
rhythm={'beta','resp'};
for rh=1:2
    
for i=1:length(SuperRat)
    myclock=tic;
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorr=SuperRat(i).trialdata.CorrIncorr10;
    
    for j=1:length(SuperRat(i).units)
        myspikes=EpochCoords(SuperRat(i).units(j).ts(:,1),snifftimes(trialCorr==1,:));
        
        for k=1:length(regions)
            % run on local beta, thats gonna be the best ref
            myLFP=SuperRat(i).([regions{k} rhythm{rh}]);
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

%                 % plot it out???                
%                 if SuperRat(i).units(j).betaRstat(k)<.01
%                     figure; ha=histogram([phasesCorr; phasesCorr+pi*2],linspace(-pi,pi*3,10+ceil(sqrt(length(phasesCorr)))),...
%                         'Normalization','probability');
%                     set(ha,'LineStyle','none','FaceColor',colors(k,:));
%                     set(gca,'XTick',[pi.*[-1:4]],'XTickLabel',{'-2\pi','','0','','2\pi'});
%                     ylabel('Probability of Spike'); xlabel('Beta Phase'); box off;
%                     title(sprintf(' %s cell in %s, %s lfp %s \n n=%.f, mean %.2f, mvl=%.2f, p=%.2e',...
%                         SuperRat(i).units(j).type, SuperRat(i).units(j).area, rhythm{rh}, regions{k}, length(phasesCorr),...
%                         SuperRat(i).units(j).betamean(k),SuperRat(i).units(j).betaMVL(k),SuperRat(i).units(j).betaRstat(k)));
%                     uiwait;
%                 end
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
for i=1:length(SuperRat)
    myclock=tic;
    snifftimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorr=SuperRat(i).trialdata.CorrIncorr10;
    
    h=waitbar(0,sprintf('Starting beta animal %s %d', SuperRat(i).name,SuperRat(i).daynum));
    for j=1:length(SuperRat(i).units)
        [~,~,~,~,myspikes]=event_spikes(SuperRat(i).units(j).ts(:,1),snifftimes(:,1),0,snifftimes(:,2)-snifftimes(:,1));
        
        for k=1:length(regions)
            myLFP=SuperRat(i).([regions{k} 'beta']);
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
%        Plot functions for both beta and resp
%
%
%
%
%


% need to run first block of claireRR analysis here
%% how many coherent cells in each region?
% now replicate claires analyses:
% local coherence and remote (cell in region 1, lfp region2)

% first, she shows a few beta locked cells....
% all odor responsive units in a given region
clear respcount betacount;
for rh=1:2
    regions={'PFC','CA1'};
    varTypes=repmat({'double'},1,6);
    allcount=table('size',[2 6],'VariableTypes',varTypes,'VariableNames',{'PFCtot','PFClocked','PFCpct','CA1tot','CA1locked','CA1pct'},'rownames',{'PFCLFP','CA1LFP'});
    %contains({SuperRat(i).units.type},'in') & ...
    
    for i=1:length(SuperRat)
        for j=1:length(regions) % cell region
            mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
                cellfun(@(a) any(a([1 2],4)<.05), {SuperRat(i).units.OdorResponsive}));
            for k=1:length(regions) % lfp region
                allcount.([regions{j} 'tot'])(k)=allcount.([regions{j} 'tot'])(k)+length(mycells);
                % add the specific lfp rstat you want here (a(k))
                allcount.([regions{j} 'locked'])(k)=allcount.([regions{j} 'locked'])(k)...
                    +sum(cellfun(@(a) a(k),{mycells.([rhythm{rh} 'Rstat'])})<.05);
            end
        end
    end
    allcount.PFCpct=allcount.PFClocked./allcount.PFCtot;
    allcount.CA1pct=allcount.CA1locked./allcount.CA1tot;
    
    eval([rhythm{rh} 'coherence=allcount;']);
    eval(['openvar ' rhythm{rh} 'coherence']);
end

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
 %% my panel c looks way different from hers

% panel E. just beta overall

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
            cellfun(@(a) any(a([1 2],4)<pcrit), {SuperRat(i).units.OdorResponsive}) & ...
            cellfun(@(a) contains(a,cctype{k}), {SuperRat(i).units.type}));       
        
            allcount.([regions{j} 'tot'])(k)=allcount.([regions{j} 'tot'])(k)+length(mycells);  % add this sess cells        
            for r=1:length(rhythm)
                allcount.([regions{j} rhythm{r}])(k)=allcount.([regions{j} rhythm{r}])(k)+sum(cellfun(@(a) a(k),{mycells.([rhythm{r} 'Rstat'])})<pcrit);
            end
        end
    end
end
 openvar('allcount')
 
 % now plot this as paired bars 
% so it will be ca1 pyr beta/resp, ca1 in beta/resp, pfc pyr beta/rsp
barpairs=[(allcount{:,[2 3]}./allcount{:,1})' , (allcount{:,[5 6]}./allcount{:,4})'];
 
hb=bar([0.9 1.1 1.4 1.6],barpairs',1,'LineStyle','none','FaceColor','flat');
for k=1:2, hb(k).CData=rhythmcolors(k,:); end
set(gca,'XTick',[0.9 1.1 1.4 1.6],'XTickLabel',{'pyr','in'});
xlabel('PFC                        CA1');
ylabel(sprintf('Proportion of \n Task Responsive Units')); box off;
legend('Beta','RR','box', 'off');
 
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
            cellfun(@(a) any(a(:,4)<pcrit),{SuperRat(i).units.OdorResponsive}));
        % for LFP location
        if ~isempty(cellpool)
            for lfploc=1:length(regions)
                % get all task responsive cells
                respp=cellfun(@(a) any(a([1 2],4)<pcrit),{cellpool.OdorResponsive});
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


pcrit=.05;
for k=1:2
    betatable=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCmean','PFCmvl','CA1mean','CA1mvl'},...
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
                    betap=cellfun(@(a) any(a([1 2],4)<pcrit),{cellpool.OdorResponsive});
                    %betap=cellfun(@(a) any(a(:,4)<pcrit),{cellpool.OdorResponsive}) &...
                    %    cell2mat(cellfun(@(a) a(betaloc)<.05, {cellpool.([rhythm{r} 'Rstat'])},'UniformOutput',false));
                    
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
            set(p,'LineStyle','none','FaceColor',colors(j,:));
            hold on;
            [rho]=circ_mean(betavals(~isnan(betavals)));
            [mvl]=circ_r(betavals(~isnan(betavals)));
            PolarArrow(rho,mvl,[],colors(i,:).*.6);
            pval=circ_rtest(betavals(~isnan(betavals)));
            
            title(sprintf('Cells %s, beta %s \n n=%d, p=%.2e',regions{j}, regions{i},length(betavals),pval));
            
        end
    end
    sgtitle(sprintf('%s',types{k}));
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
                        respp=cellfun(@(a) any(a(:,4)<pcrit),{cellpool.OdorResponsive});% &...
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
%% This is for all cells

%figure; subplot(2,1,1); histogram(bootC); subplot(2,1,2); histogram(bootI);
%linkaxes(get(gcf,'Children'),'x');

% now tabulate all the cells
% first get hpc to hpc, pfc pfc
rhythm={'beta','resp'}; clear allyy pval

saveout=0;
it=1;
for r=1:2 % rhythm
    %savefolder=uigetdir();
    allcount=table({[];[]},{[];[]},{[];[]},{[];[]},'VariableNames',{'PFCvals','PFCdprimes','CA1vals','CA1dprimes'});
    celltype='pyr';
    
    for i=1:length(SuperRat)
        for j=1:length(regions)
            mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
                contains({SuperRat(i).units.type},celltype));
            if ~isempty(mycells)
                for k=1:length(regions)
                    sessvals=cell2mat(cellfun(@(a) a(k,:), {mycells.([rhythm{r} 'PhaseDprime'])},'UniformOutput',0)');
                    allcount.([regions{j} 'vals']){k}=[allcount.([regions{j} 'vals']){k} ; sessvals(:,1:2)];
                    allcount.([regions{j} 'dprimes']){k}=[allcount.([regions{j} 'dprimes']){k}; sessvals(:,1)-sessvals(:,2)];
                end
            end
        end
    end
    
    for i=1:2
        for k=1:2
            
            figure('position',[1100,650, 300, 600 ]);  subplot(2,1,1);
            %bar([1 2],mean(allcount.([regions{i} 'vals']){k},'omitnan'),'LineStyle','none','FaceColor',colors(i,:));
            hold on;
            plot(repmat([1;2],1,length(allcount.([regions{i} 'vals']){k})),allcount.([regions{i} 'vals']){k}',...
                '-','color',rhythmcolors(r,:),'LineWidth',2);
            ylabel(sprintf('Rate Adjusted \n Vector Length'));
            box off
            xlim([0.5 2.5]); ylim([0 .6]); set(gca,'XTick',[1 2],'XTickLabel',{'Correct','Incorrect'});
            subplot(2,1,2); [y,x]=histcounts(allcount.([regions{i} 'dprimes']){k},[-.2:.02:.2],'Normalization','probability');
            barh(x(2:end)-mean(diff(x))/2,y,'LineStyle','none','BarWidth',1,'FaceColor',colors(i,:));
            hold on; plot(get(gca,'XLim'),[0 0],'k'); %repmat(mean(allcount.([regions{i} 'dprimes']){k},'omitnan'),1,2),'k');
            xlabel('Proportion of units'); ylabel('Change in MVL'); box off;
            [p]=signrank(diff(allcount.([regions{i} 'vals']){k}'));
            mymeans=mean(allcount.([regions{i} 'vals']){k},'omitnan');
            sgtitle(sprintf('All cells %s cells in %s %s LFP in %s \n n=%d,signrank=%.2e',...
                celltype, regions{i},rhythm{r}, regions{k}, length(allcount.([regions{i} 'vals']){k}'),p));
            if saveout
                figname=fullfile(savefolder,sprintf('%sCoherence allcells %s cells in %s, LFP in %s',rhythm{r},celltype,regions{i},regions{k}));
                print(figname,'-dsvg')
            end
        end
    end
    
    for i=1:2
        for k=1:2
            allyy{it}=-diff(allcount.([regions{i} 'vals']){k},1,2)';
            pval(it)=signrank(diff(allcount.([regions{i} 'vals']){k},1,2)');
            it=it+1;
        end
    end
end
% figure will be something like
% first paired points, C, IC, then a histogram beside it for difference
% so I need concatenated HPC-HPC C, IC, Dprime

%figure;
%violin(allyy);
%hold on; plot([0 9],[0 0],'k');

%% This is for all task cells, option for only locking


%figure; subplot(2,1,1); histogram(bootC); subplot(2,1,2); histogram(bootI);
%linkaxes(get(gcf,'Children'),'x');
%{
 now tabulate all the cells first get hpc to hpc, pfc pfc

changing the pcrit for the task responsiveness or the coherence metric does
not meaningfully change the main result.



%}
onlysig=0; saveout=0;
pcrit=.05;

celltype={'pyr','in'}; 
regions={'PFC','CA1'};
LFPregs={'PFC','CA1','OB'};
rhythm={'beta','resp'}; 

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
                            cellfun(@(a) any(a(1:2,4)<.05), {SuperRat(i).units.OdorResponsive}) &...
                            contains({SuperRat(i).units.type},celltype{ct}));
                    else
                        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
                            cellfun(@(a) any(a(1:2,4)<.05), {SuperRat(i).units.OdorResponsive}) &...
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
            for k=1:3 % LFP region (rows)
                try
                    figure('position',[1100,750-80*((i-1)*2+k), 300 600]);
                    subplot(2,1,1);
                    %bar([1 2],mean(allcount.([regions{i} 'vals']){k},'omitnan'),'LineStyle','none','FaceColor',colors(i,:)); hold on;
                    plot(repmat([1;2],1,length(allcount.([regions{i} 'vals']){k})),allcount.([regions{i} 'vals']){k}',...
                        '-','color',rhythmcolors(r,:),'LineWidth',3);
                    ylabel(sprintf('Rate Adjusted \n Vector Length'));
                    box off;
                    
                    xlim([0.5 2.5]); set(gca,'XTickLabel',{'Correct','Incorrect'});
                    subplot(2,1,2); [y,x]=histcounts(allcount.([regions{i} 'dprimes']){k},-.3:.05:.3,'Normalization','probability');
                    barh(x(2:end),y,'LineStyle','none','BarWidth',1,'FaceColor',colors(i,:));
                    box off; hold on; plot(get(gca,'XLim'),[0 0],'k'); %repmat(mean(allcount.([regions{i} 'dprimes']){k},'omitnan'),1,2),'k');
                    xlabel('Proportion of units'); ylabel('Change in MVL'); box off;
                    [p]=signrank(diff(allcount.([regions{i} 'vals']){k}'));
                    mymeans=mean(allcount.([regions{i} 'vals']){k},'omitnan');
                    title(sprintf('CohTaskResp %s type:%s, %s LFP in %s \n C=%.2f I=%.2f,signrank=%.2e',...
                        regions{i},celltype{ct},rhythm{r}, LFPregs{k},mymeans(1),mymeans(2),p));
                    
                   fprintf('CohTaskRe %s type:%s, %s LFP in %s \n C=%.2f I=%.2f,signrank=%.2e\n',...
                        regions{i},celltype{ct},rhythm{r}, LFPregs{k},mymeans(1),mymeans(2),p);
                    if saveout
                        % make sure you track whether these are coherent
                        % cells or not
                        figname=fullfile(savefolder,sprintf('BetaCoh Task and Coherent %s cells from %s LFP in %s',celltype,regions{i},LFPregs{k}));
                        print(figname,'-dsvg')
                    end
                end
            end
            fprintf('\n');
        end 
    end
end
set(gca,'XTick',[1 2],'YLim',[0 .6]);
%% which beta are hippocampal interneurons coherent to?  

% lets start by plot3 scatter
INcoh=[];

for i=1:length(SuperRat)
    mycells=SuperRat(i).units(contains({SuperRat(i).units.area},'CA1') &...
        cellfun(@(a) any(a(:,4)<=pcrit), {SuperRat(i).units.OdorResponsive}) &...
        cellfun(@(a) a(k)<=pcrit, {SuperRat(i).units.betaRstat}) &...
        contains({SuperRat(i).units.type},'in'));
    % lets just plot raw mvl first
    INcoh=[INcoh; cell2mat({mycells.betaMVL}')];
    
    
    
end

INcohN=INcoh-mean(INcoh,2);

%%
%
%
% Do beta coherent cells have better odor coding properties?
% Do odor coding cells have better beta coherence?
%
%

% want to aggregate all task responsive cells and then all odor selective
% cells and get their beta MVL,beta mean phase, as well as their betaphasedprime

celltable=table('Size',[2,6],'VariableTypes',repmat({'cell'},1,6),'VariableNames',...
    {'PFCbetaR','PFCbetaRho','PFCbetaDprime','CA1betaR','CA1betaRho','CA1betaDprime'},...
    'RowNames',{'responsive';'selective'});
celltype='in';

for i=1:length(SuperRat)
    for j=1:length(regions)
        % first do coders
        mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) any(a(:,4)<=pcrit), {SuperRat(i).units.OdorResponsive}) &...
            contains({SuperRat(i).units.type},celltype));
        
        celltable.([regions{j} 'betaR']){1}=[celltable.([regions{j} 'betaR']){1};...
            cellfun(@(a) a(j), {mycells.betaMVL})'];
        celltable.([regions{j} 'betaRho']){1}=[celltable.([regions{j} 'betaRho']){1};...
            cellfun(@(a) a(j), {mycells.betamean})'];
        celltable.([regions{j} 'betaDprime']){1}=[celltable.([regions{j} 'betaDprime']){1};...
            cellfun(@(a) a(j,3), {mycells.betaPhaseDprime})'];

    mycells=SuperRat(i).units(contains({SuperRat(i).units.area},regions{j}) &...
            cellfun(@(a) any(a(4)<=pcrit), {SuperRat(i).units.OdorSelective}) &...
            contains({SuperRat(i).units.type},celltype));
        
        celltable.([regions{j} 'betaR']){2}=[celltable.([regions{j} 'betaR']){2};...
            cellfun(@(a) a(j), {mycells.betaMVL})'];
        celltable.([regions{j} 'betaRho']){2}=[celltable.([regions{j} 'betaRho']){2};...
            cellfun(@(a) a(j), {mycells.betamean})'];
        celltable.([regions{j} 'betaDprime']){2}=[celltable.([regions{j} 'betaDprime']){2};...
            cellfun(@(a) a(j,3), {mycells.betaPhaseDprime})'];
    end
end


% now workup;
figure; subplot(2,1,1); histogram(celltable.PFCbetaR{1},0:.05:1);
subplot(2,1,2); histogram(celltable.PFCbetaR{2},0:.05:1);
fprintf('%.2f%% of responsive and %.2f%% of selective cells in %s are phase coherent\n',...
    nanmean(celltable.CA1betaR{1}<.05)*100,nanmean(celltable.CA1betaR{2}<.05)*100,'CA1')

fprintf('%.2f%% of responsive and %.2f%% of selective cells in %s are phase coherent\n',...
    nanmean(celltable.PFCbetaR{1}<.05)*100,nanmean(celltable.PFCbetaR{2}<.05)*100,'PFC')


figure; subplot(2,1,1);
polarhistogram(celltable.PFCbetaRho{1},15);
subplot(2,1,2);
polarhistogram(celltable.PFCbetaRho{2},15);


figure; subplot(2,1,1);
polarhistogram(celltable.CA1betaRho{1},15);
subplot(2,1,2);
polarhistogram(celltable.CA1betaRho{2},15);


figure; subplot(2,1,1); histogram(celltable.PFCbetaDprime{1},-3:.2:3);
subplot(2,1,2); histogram(celltable.PFCbetaDprime{2},-3:.2:3);












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

    