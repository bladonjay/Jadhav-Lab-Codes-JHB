%
%
%
%%


% Claire coh fig was generated using hpc tetrode with highest beta-gamma
% ratio

% claire coh fig v2 was generated using hpc and pfc tets with most units on
% them at that time


%%
%{
 % probably need a different filter here...
        % here i'm butterworthing and then hilberting.... could fiir it too
        h=waitbar(0,'gathering wavelets');
        LECspectam=abs(cwt(LECdata,bases,'cmor16-.5')); % the spectrogram for ALLLLL
        HPCspectam=abs(cwt(HPCdata,bases,'cmor16-.5'));
        
        waitbar(0,h,sprintf('Finished Wavelets, now comodulogram'));
        % then is the fast freq
        for fl=1:size(LoFreqs,1)
            % grab the lo frequencies
            [LoLphase,LoLamp]=GetLFPBand(LECdata,LFPts,LoFreqs(fl,:),1);
            [LoHphase,LoHamp]=GetLFPBand(HPCdata,LFPts,LoFreqs(fl,:),1);
            
            for fh=1:length(HiFreqs)
                
                % here you can trialwise boot by shuffling the trials and
                % correlating- if they're all the same length the
                % concatenation should line up the same
                % so to input, we'll need two transforms on each cell: 1.
                % the lo amp, and 2 the high phase that way we can only
                % pull the high phases that are also lo amp good
                loamps=cell2mat(cellfun(@(a) LoLamp(a), evinds, 'UniformOutput', false))';
                lophases=cell2mat(cellfun(@(a) LoLphase(a), evinds, 'UniformOutput', false))';
                hiamps=cell2mat(cellfun(@(a) LECspectam(fh,a), evinds, 'UniformOutput', false)');
                [~,mu,sig]=zscore(LoLamp);
                [LLcomod{fl,fh}(ses),~,~,~,~,LLnull{ses}] = CalcPhaseAmpCFC(lophases,...
                    hiamps,'LoAmps',loamps,'Cutoff',mu+sig,'runstat',100);
                
                % maybe gather the mean and std from this null, or cat all the nulls to get
                % the true dist, and then patch the real from the means
                [HHcomod{fl,fh}(ses),~,~,~,~,HHnull{ses}] = CalcPhaseAmpCFC(LoHphase(intersect(allevinds,find(OKphaseh)))',...
                    HPCspectam(fh,intersect(allevinds,find(OKphaseh)))','runstat',100);
            end
        end
%}
        
%%

% so we're first doing cross frequency coupling, not phase amp coupling


%{
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(data1,data2,params)
data1=raw first lfp
data2=raw second lfp

params includes:
params=struct('tapers',[3 5],'fs',1500,'pad',1,'Fpass',[2 40]);


%}

% the biggest choice here is which HPC lfp to use.  I am inclined to use
% the tet with the most odor responsive units or the tet with the highest
% beta to gamma power (15-30) hz to 30-45 hz ratio
%% lets see if we can get some beta to gamma ratios in tetrodes
% if there is a good spread, i'll try to relate that to place fields and
% odor fields
% per igarashi et al extended data figure 5


for i=1:length(SuperRat)
    tic
    % get behavior epochs
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10;
    odorTimes=odorTimes(diff(odorTimes')>.5,:);
    trialCorrect=trialCorrect(diff(odorTimes')>.5);
    
    % find HPC tetrodes
    myunits=find(contains({SuperRat(i).units.area},'CA1'));
    LFPtet=find(accumarray([SuperRat(i).units(myunits).tet]',1)>0);
    % generate filename here
    ratname=sprintf('%s_direct',SuperRat(i).name);
    lfpdir='EEG';
    allLFPfiles=dir(fullfile('E:\ClaireData',ratname,lfpdir));
    
    % get lfp files from today
    sessname=sprintf('%seeg%02d',SuperRat(i).name,SuperRat(i).daynum);
    todayfiles=allLFPfiles(contains({allLFPfiles.name},sessname));
    for j=1:length(todayfiles)
        todayfiles(j).tetrode=str2double(todayfiles(j).name(end-5:end-4));
    end
    % for each tetrode, accumulate all the data, get timestamps,
    % just use odor period, calculate your plots, and get ratio
    for j=1:length(LFPtet)
        % gather all files for this tetrode
        
        loadfiles=todayfiles([todayfiles.tetrode]==LFPtet(j));
        contdata=[];
        for k=1:length(loadfiles)
            lfpBit=load(fullfile(loadfiles(k).folder,loadfiles(k).name));
            tempstruct=lfpBit.eeg{SuperRat(i).daynum}{k}{LFPtet(j)};
            tempstruct.tet=LFPtet(j); tempstruct.filename=loadfiles(k).name;
            %filtered=filtereeg2(tempstruct,respfilter);
            % generate continuous data (ts, amp, instaphase, envelope)
            contdata=[contdata; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(tempstruct.data)];
        end
        
        % we'll use from 15-45, but split it in half
        params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[2 40]);
         
         eegData=[];
         for k=1:length(odorTimes)
            eegData(k,:)=contdata(contdata(:,1)>=(odorTimes(k,2)-1) & contdata(:,1)<odorTimes(k,2),2);
         end
        [S,f]=mtspectrumc(eegData',params);
        SuperRat(i).betagamma(LFPtet(j))=mean(mean(S(f<27 & f>20,:)))/mean(mean(S(f>33,:)));
        
    end
    SuperRat(i).betagamma(SuperRat(i).betagamma==0)=nan;
    fprintf('Session %d, %s %d done %d seconds \n',i, SuperRat(i).name,SuperRat(i).daynum,toc);
end


allbetagamma=cellfun(@(a) a, {SuperRat.betagamma},'UniformOutput',false);
allbetagamma2=cellfun(@(a) max(a), {SuperRat.betagamma},'UniformOutput',false);

histogram(cell2mat(allbetagamma),0:.5:6);
hold on;
histogram(cell2mat(allbetagamma2),0:.5:6);
% i think what could come next is two things:
% 1. run the cross coherence across the tetrodes, see of betagamma coupling
% pairs with beta cross coherence (may need to consider the amplitude
% confound though
% 2. cross correlate this with the number of odor coding cells on that HPC
% tetrode

%%
%
%
% Correct vs incorrect trials
%


%% overall coherence comparing correct and incorrect trials

% currently using tetrode with highest beta to gamma power in HPC

regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('DeepPink'); rgbcolormap('navy')];
combocolors=[rgbcolormap('DarkMagenta'); rgbcolormap('Chocolate'); rgbcolormap('SpringGreen')];

orders=[1 2; 1 3; 2 3];


params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[0 60]);
mindur=.75; plotIT=0;

Cperf=repmat({nan(length(SuperRat),164)},3,1); CperfX=Cperf;
Sperf=Cperf; SperfX=Cperf;

% want to do this separately for each trial I think
for i=1:length(SuperRat)
    tic
    for j=1:3
        try
            eegTime=SuperRat(i).PFCbeta(:,1);
            
            S1eeg=zscore(cell2mat({SuperRat(i).([regions{orders(j,1)} 'EEG']).data}'));
            S2eeg=zscore(cell2mat({SuperRat(i).([regions{orders(j,2)} 'EEG']).data}'));

        % now get our trial blocks
        odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
        trialCorrect=SuperRat(i).trialdata.CorrIncorr10;
        odorTimes=odorTimes(diff(odorTimes')>mindur,:);
        trialCorrect=trialCorrect(diff(odorTimes')>mindur);
        % go back from odor end
        odorTimes(:,1)=odorTimes(:,2)-mindur;
        
        C=[];S1=[]; S2=[];
        for k=1:length(odorTimes)
            okinds=eegTime>=odorTimes(k,1) & eegTime<=odorTimes(k,2);
            [C(k,:),~,~,S1(k,:),S2(k,:),f]=coherencyc(S1eeg(okinds),S2eeg(okinds),params);

        end
        %{
        if plotIT % defunct as of now
            figure; subplot(1,3,1);
            errorbar(f,nanmean(C12(trialCorrect==1,:)),SEM(C12(trialCorrect==1,:),1));
            hold on;
            errorbar(f,nanmean(C12(trialCorrect==0,:)),SEM(C12(trialCorrect==0,:),1));
            title(sprintf('HPC-PFC coherence Session %d from %s',SuperRat(i).daynum,SuperRat(i).name));
            subplot(1,3,2);
            errorbar(f,nanmean(C13(trialCorrect==1,:)),SEM(C13(trialCorrect==1,:),1));
            hold on;
            errorbar(f,nanmean(C13(trialCorrect==0,:)),SEM(C13(trialCorrect==0,:),1));
            title(sprintf('HPC-OB coherence Session %d from %s',SuperRat(i).daynum,SuperRat(i).name));
            subplot(1,3,3);
            errorbar(f,nanmean(C23(trialCorrect==1,:)),SEM(C23(trialCorrect==1,:),1));
            hold on;
            errorbar(f,nanmean(C23(trialCorrect==0,:)),SEM(C23(trialCorrect==0,:),1));
            title(sprintf('PFC-OB coherence Session %d from %s',SuperRat(i).daynum,SuperRat(i).name));
            drawnow;
        end
        %}
        % downsample here:
        minct=min([sum(trialCorrect==1) sum(trialCorrect==0)]);
        Ctrials=find(trialCorrect==1);
        Cperf{j}(i,:)=mean(C(Ctrials(randperm(length(Ctrials),minct)),:),'omitnan');
        Itrials=find(trialCorrect==0);

        CperfX{j}(i,:)=mean(C(Itrials(randperm(length(Itrials),minct)),:),'omitnan');
        
        
        Sperf{orders(j,1)}(i,:)=mean(S1(trialCorrect==1,:),'omitnan');
        SperfX{orders(j,1)}(i,:)=mean(S1(trialCorrect==0,:),'omitnan');
        Sperf{orders(j,2)}(i,:)=mean(S2(trialCorrect==1,:),'omitnan');
        SperfX{orders(j,2)}(i,:)=mean(S2(trialCorrect==0,:),'omitnan');
      
        fprintf('Session %d,  %s %d done \n',i,SuperRat(i).name,SuperRat(i).daynum);
    catch
        fprintf('Session %d with %s didnt work \n',SuperRat(i).daynum,SuperRat(i).name);
        end
    end
end

% Tried this with HPC tetrode with most cells vs PFC tet with most cells
% results did not look good

% then tried this with HPC tetrode with highest beta-gamma ratio with PFC
% tet with most cells, that did not look good either

% I think we need more tets with distal CA1, I have a feeling most of these
% recordings do not have distal tets and therefore are not useful for this
% analysis

%% now plotting functions


% thinking a full spectrogram with patch, then inset is...
% a histogram of shaded regions differences by session, to show that across
% sessions the effect is consistent (or isnt!)

% first coherence beta is 20-30, rr is 7-8 with no movement


figure;
clear pl;
for k=1:3
    % start with patches in back
    ylims=[-.2 .2];
    % change range colors to our colors
    if k==1
        patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
        hold on; patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
    end
    pl(k)=plot(f,mean(Cperf{k}-CperfX{k},'omitnan'),'color',combocolors(k,:),'LineWidth',2);
    patch([f fliplr(f)]', [mean(Cperf{k}-CperfX{k},'omitnan')+SEM(Cperf{k}-CperfX{k},1)...
        fliplr(mean(Cperf{k}-CperfX{k},'omitnan')-SEM(Cperf{k}-CperfX{k},1))]',...
        combocolors(k,:),'LineStyle','none','FaceAlpha',.5);
    
    ylim(ylims);
    xlabel('Frequency, Hz');
    ylabel(sprintf('Correct-Incorrect \n Cross Coherence'));
end

plot([0 60],[0 0],'k','LineWidth',3);
legend(pl,{'CA1-PFC','CA1-OB','PFC-OB'});


%% the old way, separate plots

for i=1:3
    figure;
    % inset
    subplot(1,2,1);
    % get mean rr amp across
    RRmeans=mean(Cperf{i}(:,f>=7 & f<=8),2,'omitnan');
    RRPremeans=mean(CperfX{i}(:,f>=7 & f<=8),2,'omitnan');
    [a,b]=histcounts(RRmeans-RRPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(1,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(RRmeans-RRPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(1,:),...
        'LineWidth',2); box off; ylabel('Counts'); xlabel('Change in CA1-PFC RR Coherence')
    title(sprintf('p=%.2e',signrank(RRmeans-RRPremeans)));
    
    
    subplot(1,2,2);
    % get mean beta amp across
    betameans=mean(Cperf{i}(:,f>=20 & f<=30),2,'omitnan');
    betaPremeans=mean(CperfX{i}(:,f>=20 & f<=30),2,'omitnan');
    [a,b]=histcounts(betameans-betaPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(2,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(betameans-betaPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(2,:),...
        'LineWidth',2); box off; ylabel('Counts'); xlabel('Change CA1-PFC in Beta Coherence')
    title(sprintf('p=%.2e',signrank(betameans-betaPremeans)));
    sgtitle(sprintf('%s %s coherence',regions{orders(i,1)},regions{orders(i,2)}));
end


%% Now Power results (have to multiply by freq

figure;
clear pl;
for k=1:3
    % start with patches in back
    ylims=[-.15 .15];
    % change range colors to our colors
    if k==1
        patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
        hold on; patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
    end
    pl(k)=plot(f,mean(Sperf{k}-SperfX{k},'omitnan').*f,'color',colors(k,:),'LineWidth',2);
    patch([f fliplr(f)]', [mean(Sperf{k}-SperfX{k},'omitnan').*f+SEM(Sperf{k}-SperfX{k},1).*f...
        fliplr(mean(Sperf{k}-SperfX{k},'omitnan').*f-SEM(Sperf{k}-SperfX{k},1).*f)]',...
        colors(k,:),'LineStyle','none','FaceAlpha',.5);
    
    ylim(ylims);
    xlabel('Frequency, Hz');
    ylabel(sprintf('Correct-Incorrect \n Power Change (db)'));
end
plot([0 60],[0 0],'k','LineWidth',3);
legend(pl,regions);
% and the histograms if you want them...


%%
for k=1:3
    figure;
    % inset
    subplot(1,2,1);
    % get mean rr amp across
    RRmeans=mean(Sperf{k}(:,f>=7 & f<=8),2,'omitnan');
    RRPremeans=mean(SperfX{k}(:,f>=7 & f<=8),2,'omitnan');
    [a,b]=histcounts(RRmeans-RRPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(1,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(RRmeans-RRPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(1,:),...
        'LineWidth',2); box off; ylabel('Counts');
    xlabel(sprintf('Change in %s RR Power',regions{k}));
    title(sprintf('p=%.2e',signrank(RRmeans-RRPremeans)));
    
    
    subplot(1,2,2);
    % get mean beta amp across
    betameans=mean(Sperf{k}(:,f>=20 & f<=30),2,'omitnan');
    betaPremeans=mean(SperfX{k}(:,f>=20 & f<=30),2,'omitnan');
    [a,b]=histcounts(betameans-betaPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(2,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(betameans-betaPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(2,:),...
        'LineWidth',2); box off; ylabel('Counts');
    xlabel(sprintf('Change in %s Beta Power (db)',regions{k}));
    title(sprintf('p=%.2e',signrank(betameans-betaPremeans)));
end

%%

%
%
% pre odor vs odor period
%
%
%



%% lets now try power following onset (or before offset)

regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')];
orders=[1 2; 1 3; 2 3];

params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[0 60]);
mindur=.75;
Csess=repmat({nan(length(SuperRat),164)},3,1); % this is based on the mindur and the fpass
preCsess=Csess;
orders=[1 2; 1 3; 2 3]; %
Ssess=Csess; preSsess=Csess;



% want to do this separately for each trial I think
for i=1:length(SuperRat)
    tic
    for j=1:3
        try
            eegTime=SuperRat(i).PFCbeta(:,1);
            
            S1eeg=zscore(cell2mat({SuperRat(i).([regions{orders(j,1)} 'EEG']).data}'));
            S2eeg=zscore(cell2mat({SuperRat(i).([regions{orders(j,2)} 'EEG']).data}'));
            % now get our trial blocks
            
            odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend]; % full poke time
            odorTimes=odorTimes(diff(odorTimes')>mindur,:);
            trialCorrect=SuperRat(i).trialdata.CorrIncorr10(diff(odorTimes')>mindur);
            
            odorTimes2=[odorTimes(:,1)-mindur odorTimes(:,1)]; % backwards from odor start
            odorTimes(:,1)=odorTimes(:,2)-mindur; % go backwards from odor end
            clear C S1 S2 preC preS1 preS2;
            for k=1:length(odorTimes)
                okinds=eegTime>=odorTimes(k,1) & eegTime<=odorTimes(k,2);
                [C(k,:),~,~,S1(k,:),S2(k,:),f]=coherencyc(S1eeg(okinds),S2eeg(okinds),params);
                okinds=eegTime>=odorTimes2(k,1) & eegTime<=odorTimes2(k,2);
                [preC(k,:),~,~,preS1(k,:),preS2(k,:),f]=coherencyc(S1eeg(okinds),S2eeg(okinds),params);
                
            end
            % plot each session out
            %{
    figure;
    subplot(1,3,1);
    errorbar(f,nanmean(S1(trialCorrect==1,:)),SEM(S1(trialCorrect==1,:),1));
    hold on;
    errorbar(f,nanmean(preS1(trialCorrect==1,:)),SEM(preS1(trialCorrect==1,:),1));
    title('CA1 Power'); legend('poke','pre-poke');
    %plot(f,zscore(C(trialCorrect~=1,:),1,2),'r');
    subplot(1,3,2);
    errorbar(f,nanmean(S2(trialCorrect==1,:)),SEM(S2(trialCorrect==1,:),1));
    hold on;
    errorbar(f,nanmean(preS2(trialCorrect==1,:)),SEM(preS2(trialCorrect==1,:),1));
    title('PFC Power'); legend('poke','pre-poke');
    subplot(1,3,3);
    errorbar(f,nanmean(C(trialCorrect==1,:)),SEM(C(trialCorrect==1,:),1));
    hold on;
    errorbar(f,nanmean(preC(trialCorrect==1,:)),SEM(preC(trialCorrect==1,:),1));
    title('PFC-CA1 Coherence'); legend('poke','pre-poke');
    drawnow;
            %}
            
            % or save the data into a sessionwise struct
            Csess{j}(i,:)=mean(C); preCsess{j}(i,:)=mean(preC);
            Ssess{orders(j,1)}(i,:)=mean(S1); preSsess{orders(j,1)}(i,:)=mean(preS1);
            Ssess{orders(j,2)}(i,:)=mean(S2); preSsess{orders(j,2)}(i,:)=mean(preS2);
            
            fprintf('Session %d %s %d took %d seconds \n',i,SuperRat(i).name,SuperRat(i).daynum,round(toc))
            
        catch
            fprintf('Session %d, %d with %s didnt work \n',i, SuperRat(i).daynum,SuperRat(i).name)
        end
    end
end

%% now plot Coherence Results

% plot before and after reward onset here

% thinking a full spectrogram with patch, then inset is...
% a histogram of shaded regions differences by session, to show that across
% sessions the effect is consistent (or isnt!)

% first coherence beta is 20-30, rr is 7-8 with no movement
regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
combocolors=[rgbcolormap('DarkMagenta'); rgbcolormap('Chocolate'); rgbcolormap('SpringGreen')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')];
figure('position',[1000 910 910 450]);  
for i=1:3
    % start with patches in back
    ylims=[-.1 .3];
    % change range colors to our colors
    if i==1
        patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
        hold on; patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
        
        p=plot(f,mean(Csess{i}-preCsess{i},'omitnan'),'Color',combocolors(i,:),'LineWidth',2);
        patch([f fliplr(f)]', [mean(Csess{i}-preCsess{i},'omitnan')+SEM(Csess{i}-preCsess{i},1)...
            fliplr(mean(Csess{i}-preCsess{i},'omitnan')-SEM(Csess{i}-preCsess{i},1))]',...
            combocolors(i,:),'LineStyle','none','FaceAlpha',.5);
    else
        p(i)=plot(f,mean(Csess{i}-preCsess{i},'omitnan'),'Color',combocolors(i,:),'LineWidth',2);
        patch([f fliplr(f)]', [mean(Csess{i}-preCsess{i},'omitnan')+SEM(Csess{i}-preCsess{i},1)...
            fliplr(mean(Csess{i}-preCsess{i},'omitnan')-SEM(Csess{i}-preCsess{i},1))]',...
            combocolors(i,:),'LineStyle','none','FaceAlpha',.5);
        ylim(ylims);
    end
    xlabel('Frequency, Hz');
    ylabel(sprintf('Pre-to-Post Odor Onset\nChange in Coherence'));
end
plot([0 60],[0 0],'k','LineWidth',3);
legend(p,{'CA1-PFC','CA1-OB','PFC-OB'});

%% and do separate plots here
for i=1:3
figure;
% inset
subplot(1,2,1);
% get mean rr amp across
RRmeans=mean(Csess{i}(:,f>=7 & f<=8),2,'omitnan');
RRPremeans=mean(preCsess{i}(:,f>=7 & f<=8),2,'omitnan');
[a,b]=histcounts(RRmeans-RRPremeans,10);
bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
    'FaceColor',rhythmcolors(1,:),'FaceAlpha',.7);
hold on; plot(repmat(mean(RRmeans-RRPremeans,'omitnan'),1,2),[0 7],'Color',rhythmcolors(1,:),...
    'LineWidth',2); box off; ylabel('Counts'); xlabel('Change in RR Coherence')
title(sprintf('signrank p=%.2e',signrank(RRmeans-RRPremeans)));


subplot(1,2,2);
% get mean beta amp across
betameans=mean(Csess{i}(:,f>=20 & f<=30),2,'omitnan');
betaPremeans=mean(preCsess{i}(:,f>=20 & f<=30),2,'omitnan');
[a,b]=histcounts(betameans-betaPremeans,10);
bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
    'FaceColor',rhythmcolors(2,:),'FaceAlpha',.7);
hold on; plot(repmat(mean(betameans-betaPremeans,'omitnan'),1,2),[0 9],'Color',rhythmcolors(2,:),...
    'LineWidth',2); box off; ylabel('Counts'); xlabel('Change in Beta Coherence')
title(sprintf('p=%.2e',signrank(betameans-betaPremeans)));
sgtitle(sprintf('%s -%s coherence',regions{orders(i,1)},regions{orders(i,2)}));
end


%%
% plot difference in coherence between before and after for all three
% regions
%{
for i=1:3
  figure; 
subplot(2,2,[1 2]);
% start with patches in back
ylims=[.4 .7];
% change range colors to our colors
patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
hold on;
patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);

p=plot(f,nanmean(Csess{i}-preCsess{i}),'b','LineWidth',2);
patch([f fliplr(f)]', [mean(Csess{i},'omitnan')+SEM(Csess{i},1) fliplr(mean(Csess{i},'omitnan')-SEM(Csess{i},1))]',...
    'b','LineStyle','none','FaceAlpha',.5);
p(2)=plot(f,nanmean(preCsess{i}),'k','LineWidth',2);
patch([f fliplr(f)]', [mean(preCsess{i},'omitnan')+SEM(preCsess{i},1) fliplr(mean(preCsess{i},'omitnan')-SEM(preCsess{i},1))]',...
    'k','LineStyle','none','FaceAlpha',.5);
ylim(ylims);
legend(p,'Odor Period','Pre-odor'); xlabel('Frequency, Hz');
ylabel(sprintf('%s-%s Coherence',regions{orders(i,1)},regions{orders(i,2)}));



% inset
subplot(2,2,3);
% get mean rr amp across
RRmeans=mean(Csess{i}(:,f>=7 & f<=8),2,'omitnan');
RRPremeans=mean(preCsess{i}(:,f>=7 & f<=8),2,'omitnan');
[a,b]=histcounts(RRmeans-RRPremeans,10);
bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
    'FaceColor',rhythmcolors(1,:),'FaceAlpha',.7);
hold on; plot(repmat(mean(RRmeans-RRPremeans,'omitnan'),1,2),[0 7],'Color',rhythmcolors(1,:),...
    'LineWidth',2); box off; ylabel('Counts'); xlabel('Change in RR Coherence')
title(sprintf('p=%.2e',signrank(RRmeans-RRPremeans)));


subplot(2,2,4);
% get mean beta amp across
betameans=mean(Csess{i}(:,f>=20 & f<=30),2,'omitnan');
betaPremeans=mean(preCsess{i}(:,f>=20 & f<=30),2,'omitnan');
[a,b]=histcounts(betameans-betaPremeans,10);
bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
    'FaceColor',rhythmcolors(2,:),'FaceAlpha',.7);
hold on; plot(repmat(mean(betameans-betaPremeans,'omitnan'),1,2),[0 9],'Color',rhythmcolors(2,:),...
    'LineWidth',2); box off; ylabel('Counts'); xlabel('Change in Beta Coherence')
title(sprintf('p=%.2e',signrank(betameans-betaPremeans)));

end
%}
%% Now Power results

% all three on same plot
figure('position',[1000 910 910 450]);
clear pl;
for k=1:3
    S1sessdb=Ssess{k}.*f; preS1sessdb=preSsess{k}.*f;

    % start with patches in back
    ylims=[-.15 .45];
    if k==1
        % change range colors to our colors
        patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
        hold on; patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
    end
    pl(k)=plot(f,mean(S1sessdb-preS1sessdb,'omitnan'),'Color',colors(k,:),'LineWidth',2);
    patch([f fliplr(f)]', [mean(S1sessdb-preS1sessdb,'omitnan')+SEM(S1sessdb-preS1sessdb,1)...
        fliplr(mean(S1sessdb-preS1sessdb,'omitnan')-SEM(S1sessdb-preS1sessdb,1))]',...
        colors(k,:),'LineStyle','none','FaceAlpha',.5);
   
    ylim(ylims);
    xlabel('Frequency, Hz');
    ylabel(sprintf('Pre-to-Post Odor Onset \n Change in Power (db)'));
end
plot([0 60],[0 0],'k','LineWidth',3);
legend(pl,regions);

%% on separate plots
maxpow=[.35 .35 .6];
for k=1:3
    S1sessdb=Ssess{k}.*f; preS1sessdb=preSsess{k}.*f;
    figure;
    subplot(2,2,1:2);
    % start with patches in back
    ylims=[0 maxpow(k)];
    % change range colors to our colors
    patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
    hold on;
    patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
    p=plot(f,nanmean(S1sessdb),'Color',colors(k,:),'LineWidth',2);
    patch([f fliplr(f)]', [mean(S1sessdb,'omitnan')+SEM(S1sessdb,1) fliplr(mean(S1sessdb,'omitnan')-SEM(S1sessdb,1))]',...
        colors(k,:),'LineStyle','none','FaceAlpha',.5);
    p(2)=plot(f,nanmean(preS1sessdb),'k','LineWidth',2);
    patch([f fliplr(f)]', [mean(preS1sessdb,'omitnan')+SEM(preS1sessdb,1) fliplr(mean(preS1sessdb,'omitnan')-SEM(preS1sessdb,1))]',...
        'k','LineStyle','none','FaceAlpha',.5);
    ylim(ylims);
    legend(p,'Odor Period','Pre-odor'); xlabel('Frequency, Hz');
    ylabel(sprintf('%s Power (db)',regions{k}));
    
    
    
    
    % inset
    subplot(2,2,3);
    % get mean rr amp across
    RRmeans=mean(Ssess{k}(:,f>=7 & f<=8),2,'omitnan');
    RRPremeans=mean(preSsess{k}(:,f>=7 & f<=8),2,'omitnan');
    [a,b]=histcounts(RRmeans-RRPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(1,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(RRmeans-RRPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(1,:),...
        'LineWidth',2); box off; ylabel('Counts');
    xlabel(sprintf('Change in %s RR Power',regions{k}))
    title(sprintf('p=%.2e',signrank(RRmeans-RRPremeans)));
    
    
    subplot(2,2,4);
    % get mean beta amp across
    betameans=mean(Ssess{k}(:,f>=20 & f<=30),2,'omitnan');
    betaPremeans=mean(preSsess{k}(:,f>=20 & f<=30),2,'omitnan');
    [a,b]=histcounts(betameans-betaPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(2,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(betameans-betaPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(2,:),...
        'LineWidth',2); box off; ylabel('Counts');
    xlabel(sprintf('Change in %s Beta Power',regions{k}));
    title(sprintf('p=%.2e',signrank(betameans-betaPremeans)));
end

%% Now try odor period to reward period


%
%
% odor Period to Reward Period
%
%
%



%% match to seconds preceding release

params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[0 60]);
mindur=.75;

Csess=repmat({nan(length(SuperRat),164)},3,1); % this is based on the mindur and the fpass
CsessR=Csess;
Ssess=Csess; SsessR=Csess;

orders=[1 2; 1 3; 2 3]; %


% want to do this separately for each trial I think
for i=1:length(SuperRat)
    tic
    for j=1:3
        try
            eegTime=SuperRat(i).PFCbeta(:,1);
            
            S1eeg=zscore(cell2mat({SuperRat(i).([regions{orders(j,1)} 'EEG']).data}'));
            S2eeg=zscore(cell2mat({SuperRat(i).([regions{orders(j,2)} 'EEG']).data}'));
            % now get our trial blocks
            
            odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend]; % full poke time
            if isfield(SuperRat(i).trialdata,'rewardend')
                rewardTimes=[SuperRat(i).trialdata.rewardstart SuperRat(i).trialdata.rewardend]; % full poke time
            elseif isfield(SuperRat(i).trialdata,'rewarden')
                rewardTimes=[SuperRat(i).trialdata.rewardstart SuperRat(i).trialdata.rewarden]; % full poke time
            end
            oktrials=SuperRat(i).trialdata.CorrIncorr10' & diff(odorTimes')>mindur;
            odorTimes=odorTimes(oktrials,:); rewardTimes=rewardTimes(oktrials,:);
            odorTimes(:,1)=odorTimes(:,2)-mindur; % go backwards from odor end
            rewardTimes(:,1)=rewardTimes(:,2)-mindur;
            clear C S1 S2 preC preS1 preS2;
            for k=1:length(odorTimes)
                okinds=eegTime>=odorTimes(k,1) & eegTime<=odorTimes(k,2);
                [C(k,:),~,~,S1(k,:),S2(k,:),f]=coherencyc(S1eeg(okinds),S2eeg(okinds),params);
                okinds=eegTime>=rewardTimes(k,1) & eegTime<=rewardTimes(k,2);
                [preC(k,:),~,~,preS1(k,:),preS2(k,:),f]=coherencyc(S1eeg(okinds),S2eeg(okinds),params);
                
            end

            
            % or save the data into a sessionwise struct
            Csess{j}(i,:)=mean(C); CsessR{j}(i,:)=mean(preC);
            Ssess{orders(j,1)}(i,:)=mean(S1); SsessR{orders(j,1)}(i,:)=mean(preS1);
            Ssess{orders(j,2)}(i,:)=mean(S2); SsessR{orders(j,2)}(i,:)=mean(preS2);
            
            fprintf('Session %d %s %d took %d seconds \n',i,SuperRat(i).name,SuperRat(i).daynum,round(toc))
            
        catch
            fprintf('Session %d, %d with %s didnt work \n',i, SuperRat(i).daynum,SuperRat(i).name)
        end
    end
end
%% power odor period to reard period

% all three on same plot
figure('position',[1000 910 910 450]);
clear pl;
for k=1:3
    S1sessdb=Ssess{k}.*f; preS1sessdb=SsessR{k}.*f;

    % start with patches in back
    ylims=[-.15 .65];
    if k==1
        % change range colors to our colors
        patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
        hold on; patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
    end
    pl(k)=plot(f,mean(S1sessdb-preS1sessdb,'omitnan'),'Color',colors(k,:),'LineWidth',2);
    patch([f fliplr(f)]', [mean(S1sessdb-preS1sessdb,'omitnan')+SEM(S1sessdb-preS1sessdb,1)...
        fliplr(mean(S1sessdb-preS1sessdb,'omitnan')-SEM(S1sessdb-preS1sessdb,1))]',...
        colors(k,:),'LineStyle','none','FaceAlpha',.5);
   
    ylim(ylims);
    xlabel('Frequency, Hz');
    ylabel(sprintf('Odor-Reward \n Change in Power (db)'));
end
plot([0 60],[0 0],'k','LineWidth',3);
legend(pl,regions);

%% on separate plots
maxpow=[.35 .35 .6];
for k=1:3
    S1sessdb=Ssess{k}.*f; preS1sessdb=SsessR{k}.*f;
    figure;
    subplot(2,2,1:2);
    % start with patches in back
    ylims=[0 maxpow(k)];
    % change range colors to our colors
    patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
    hold on;
    patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
    p=plot(f,nanmean(S1sessdb),'Color',colors(k,:),'LineWidth',2);
    patch([f fliplr(f)]', [mean(S1sessdb,'omitnan')+SEM(S1sessdb,1) fliplr(mean(S1sessdb,'omitnan')-SEM(S1sessdb,1))]',...
        colors(k,:),'LineStyle','none','FaceAlpha',.5);
    p(2)=plot(f,nanmean(preS1sessdb),'k','LineWidth',2);
    patch([f fliplr(f)]', [mean(preS1sessdb,'omitnan')+SEM(preS1sessdb,1) fliplr(mean(preS1sessdb,'omitnan')-SEM(preS1sessdb,1))]',...
        'k','LineStyle','none','FaceAlpha',.5);
    ylim(ylims);
    legend(p,'Odor Period','Pre-odor'); xlabel('Frequency, Hz');
    ylabel(sprintf('%s Power (db)',regions{k}));
    
    
    
    
    % inset
    subplot(2,2,3);
    % get mean rr amp across
    RRmeans=mean(Ssess{k}(:,f>=7 & f<=8),2,'omitnan');
    RRPremeans=mean(SsessR{k}(:,f>=7 & f<=8),2,'omitnan');
    [a,b]=histcounts(RRmeans-RRPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(1,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(RRmeans-RRPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(1,:),...
        'LineWidth',2); box off; ylabel('Counts');
    xlabel(sprintf('Change in %s RR Power',regions{k}))
    title(sprintf('p=%.2e',signrank(RRmeans-RRPremeans)));
    
    
    subplot(2,2,4);
    % get mean beta amp across
    betameans=mean(Ssess{k}(:,f>=20 & f<=30),2,'omitnan');
    betaPremeans=mean(SsessR{k}(:,f>=20 & f<=30),2,'omitnan');
    [a,b]=histcounts(betameans-betaPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(2,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(betameans-betaPremeans,'omitnan'),1,2),[0 max(a)],'Color',rhythmcolors(2,:),...
        'LineWidth',2); box off; ylabel('Counts');
    xlabel(sprintf('Change in %s Beta Power',regions{k}));
    title(sprintf('p=%.2e',signrank(betameans-betaPremeans)));
end

%% coherence

% plot before and after reward onset here

% thinking a full spectrogram with patch, then inset is...
% a histogram of shaded regions differences by session, to show that across
% sessions the effect is consistent (or isnt!)

% first coherence beta is 20-30, rr is 7-8 with no movement
regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
combocolors=[rgbcolormap('DarkMagenta'); rgbcolormap('Chocolate'); rgbcolormap('SpringGreen')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')];
figure('position',[1000 910 910 450]);  
for i=1:3
    % start with patches in back
    ylims=[-.1 .3];
    % change range colors to our colors
    if i==1
        patch([7 8 8 7],[linearize(repmat(ylims,2,1))],rhythmcolors(1,:),'LineStyle','none','FaceAlpha',.3);
        hold on; patch([20 30 30 20],[linearize(repmat(ylims,2,1))],rhythmcolors(2,:),'LineStyle','none','FaceAlpha',.3);
        
        p=plot(f,mean(Csess{i}-CsessR{i},'omitnan'),'Color',combocolors(i,:),'LineWidth',2);
        patch([f fliplr(f)]', [mean(Csess{i}-CsessR{i},'omitnan')+SEM(Csess{i}-CsessR{i},1)...
            fliplr(mean(Csess{i}-CsessR{i},'omitnan')-SEM(Csess{i}-CsessR{i},1))]',...
            combocolors(i,:),'LineStyle','none','FaceAlpha',.5);
    else
        p(i)=plot(f,mean(Csess{i}-CsessR{i},'omitnan'),'Color',combocolors(i,:),'LineWidth',2);
        patch([f fliplr(f)]', [mean(Csess{i}-CsessR{i},'omitnan')+SEM(Csess{i}-CsessR{i},1)...
            fliplr(mean(Csess{i}-CsessR{i},'omitnan')-SEM(Csess{i}-CsessR{i},1))]',...
            combocolors(i,:),'LineStyle','none','FaceAlpha',.5);
        ylim(ylims);
    end
    xlabel('Frequency, Hz');
    ylabel(sprintf('Odor-Reward\nChange in Coherence'));
end
plot([0 60],[0 0],'k','LineWidth',3);
legend(p,{'CA1-PFC','CA1-OB','PFC-OB'});



%% separate plots and statistiscs
for i=1:3
    figure;
    % inset
    subplot(1,2,1);
    % get mean rr amp across
    RRmeans=mean(Csess{i}(:,f>=7 & f<=8),2,'omitnan');
    RRPremeans=mean(CsessR{i}(:,f>=7 & f<=8),2,'omitnan');
    [a,b]=histcounts(RRmeans-RRPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(1,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(RRmeans-RRPremeans,'omitnan'),1,2),[0 7],'Color',rhythmcolors(1,:),...
        'LineWidth',2); box off; ylabel('Counts'); xlabel('Change in RR Coherence')
    title(sprintf('signrank p=%.2e',signrank(RRmeans-RRPremeans)));
    
    
    subplot(1,2,2);
    % get mean beta amp across
    betameans=mean(Csess{i}(:,f>=20 & f<=30),2,'omitnan');
    betaPremeans=mean(CsessR{i}(:,f>=20 & f<=30),2,'omitnan');
    [a,b]=histcounts(betameans-betaPremeans,10);
    bar(mean([b(1:end-1); b(2:end)]),a,1,'LineStyle','none',...
        'FaceColor',rhythmcolors(2,:),'FaceAlpha',.7);
    hold on; plot(repmat(mean(betameans-betaPremeans,'omitnan'),1,2),[0 9],'Color',rhythmcolors(2,:),...
        'LineWidth',2); box off; ylabel('Counts'); xlabel('Change in Beta Coherence')
    title(sprintf('p=%.2e',signrank(betameans-betaPremeans)));
    sgtitle(sprintf('%s -%s coherence odor to reward',regions{orders(i,1)},regions{orders(i,2)}));
end


%% now calculate phase amp CFC.


%{
igarashi found significant rr-beta phase amplitude coupling during odor
sampling, I have not ofund that to be true. I actually think there is a
shift from rr to beta as the animal preps his action.
  
I'm not sure why they find this coherence now...

So the question is
why is the beta coherence not as good for incorrect trials if that is
likely the evidence of the post-decision period?  


%}

%
%
%
%
%
%
%
%
%%


% crucial function is calcPhaseAmpCFC

fnq=1500/2; % nyquist
wname='cmor16-.5'; % a 16 x envelope of a morlet at 0.5 hz central frequency
freqs=logspace(.5,2.3,120); % log spaced vector to ~250 hz
%freqs=linspace(2,150,80);
bases=fnq./freqs; % 500 is fnq, so this is base FNQ
% I think i can run this trialwise.... but i may need to do whole
% session...

% in hz

phaserange=[0 40];
phaseInds=find(freqs>phaserange(1) & freqs<phaserange(2));
phaseF=freqs(phaseInds);
amprange=[12 200];
ampInds=find(freqs>amprange(1) & freqs<amprange(2));
ampF=freqs(ampInds);
% start with one session, do pre-post
inspect=repmat({},length(SuperRat),3);
xspect=repmat({},length(SuperRat),3);

for i=1:length(SuperRat)
    eegTime=SuperRat(i).PFCbeta(:,1);
    
    % now get our trial blocks
    mindur=0.5; maxdur=2; buffer=.5;
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend]; % full poke time
    odorTimes=odorTimes(diff(odorTimes')>mindur & diff(odorTimes')<maxdur,:);
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10(diff(odorTimes')>mindur & diff(odorTimes')<maxdur);
    
    odorTimes2=[odorTimes(:,1)-mindur odorTimes(:,1)]; % backwards from odor start
    odorTimes(:,1)=odorTimes(:,2)-mindur; % go backwards from odor end
    
    % so maybe run this trialwise, but add a buffer, and then remove the buffer
    % afterwards.  The buffer should be something like twice the duration of
    % the slowest oscillation i want to see... so 5hz.. 2.5 hz, or 4/10 of a
    % second, say 500 msec round up
    
    % this is probably too short a window... need to get a bigger one and
    % compress...
    [~,~,~,inds,~,times]=event_spikes(eegTime,odorTimes(:,1),1,2);
    %mindur=min(cellfun(@(a) length(a), inds));
    % now we have the raw data, and we can perform our cwts
    
    rawEEG=repmat({},3,1);
    for k=1:3
        rawEEG{k}=zscore(cell2mat({SuperRat(i).([regions{k} 'EEG']).data}'));
    end
    
    if plotIT==1, figure; end
    
    for k=2 %1:length(regions) % first set of spects (same)
        
        phaseDist=nan(length(apr),50);
        PAcomod=nan(length(freqs),length(freqs));
        PAcomod2=PAcomod;
        trspectam={};
        for j=1:length(inds)
            % for each trial, get a full cross spectrogram, save that spect
            temp=cwt(rawEEG{k}(inds{j}),bases,wname); % the eeg for this trial
            trspectam{j}=temp(:,times{j}>=0 & times{j}<=.5); % just use inside vals
        end
        minsz=min(cellfun(@(a) size(a,2), trspectam));
        trspectam=cellfun(@(a) a(:,1:minsz), trspectam,'UniformOutput',false);
        spectam=cell2mat(trspectam(trialCorrect'==1));
        spectamX=cell2mat(trspectam(trialCorrect'==0));

        RRphasedist=nan(length(ampInds),50);
        % now cfc that and asave it
        for phaseidx = 1:length(phaseInds) % low freq
            for apr = ampInds(ampInds>phaseidx) % high freq
                [PAcomod(apr,phaseidx)] = CalcPhaseAmpCFC(angle(spectam(phaseidx,:)),...
                    abs(spectam(apr,:)));
                [PAcomod2(apr,phaseidx)] = CalcPhaseAmpCFC(angle(spectamX(phaseidx,:)),...
                    abs(spectamX(apr,:)));
                
                %if phaseidx==53
                %    [PAcomod(apr,phaseidx),RRphasedist(apr,:)] = CalcPhaseAmpCFC(angle(spectam(phaseidx,:)),...
                %        abs(spectam(apr,:)));
                %end
            end
        end
        inspect{i,k}=PAcomod;
        Xspect{i,k}=PAcomod2;
        if plotIT==1
            subplot(3,2,(k*2)-1);
            imagesc(freqs,freqs,PAcomod)
            set(gca,'Ydir','normal','YScale','log','XSCale','log');
            ylim(amprange); xlim(phaserange);
            title(sprintf('%s phaseampcomod',regions{k}));
            xlabel('Phase freq'); ylabel('amp freq');
            subplot(3,2,k*2);
            imagesc(linspace(-pi,pi,50),ampF,RRphasedist(ampInds(1):end,:))
            set(gca,'Ydir','normal','YScale','log')
            xlabel('Phase of 7Hz'); ylabel('amp freq');
        end
        fprintf('%d session %s %d done \n',i,SuperRat(i).name,SuperRat(i).daynum);
    end
end
%% now plot the data in the aggregate
figure;
for i=1:3
    subplot(3,1,i);
    newcell=cell2mat(permute(inspect(:,i),[3 2 1]));
    imagesc(freqs,freqs,mean(newcell,3,'omitnan'));
    set(gca,'Ydir','normal','YScale','log','XSCale','log');
    ylim(amprange); xlim(phaserange);
    set(gca,'XTick',[2 4 8 12 16 20 30 40],'YTick',[12 15 20 30 50 80 120 160 200]);
    title([regions{i}]);
end

% plot the difference between correct and incorrect in the aggregate
diffscores=cellfun(@(a,b) (a-b), inspect(:,2),Xspect(:,2),'UniformOutput',false);

figure;
i=2;

newcell=cell2mat(permute(diffscores,[3 2 1]));
imagesc(freqs,freqs,mean(newcell,3,'omitnan'));
set(gca,'Ydir','normal','YScale','log','XSCale','log');
ylim(amprange); xlim(phaserange);
set(gca,'XTick',[2 4 8 12 16 20 30 40],'YTick',[12 15 20 30 50 80 120 160 200]);
title([regions{i}]);

%%
% inspect a few trials:
for j=8:2:22
    
    figure; subplot(2,1,1);
    plot(times{j},rawEEG{k}(inds{j}))
    subplot(2,1,2); imagesc(times{j},freqs,abs(trspectam{j}))
    set(gca,'ydir','normal','yscale','log','YTick',[2 5 8 12 20 30 50 80])
    ylim([5 80]);
end



LECspectam=cwt(zscore(LECdata),bases,wname); % the spectrogram for ALLLLL

%%
%
%
%
% replicating claires beta coherence onset analysis
% the reviewers didnt like that 
%
%
%
%
%
%
%
