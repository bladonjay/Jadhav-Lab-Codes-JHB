%%
%
%
%
%
% Calculating phase offsets for beta and RR during odor sampling
%
%
%
%
%

% observation block:
%{
- first observation is that locked to sniff end you get a small peak at end
in beta offset in the positive direction
- if you lock to sniff start, you get a positive swing then a massive
negative swing
- There are definitiely fluctuations locked to sniff start and sniff end,
but its hard to suss out whats going on in aggregate.

%}


regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')]; % BETA IS BLUE
types={'pyr','in'};
rhythm={'beta','resp'};


%%
%
%   observe single sessions here
%
%

% i think the gist of this is that the beta power just keeps going up the
% longer the animal is in the odor port...
pairs=[1 2; 1 3; 2 3];
span=[-.7 .7]; % secbefore to secafter
dur=1500*diff(span);
for i=[4 15 22] %length(SuperRat)
    % get odor starts
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    odorTimes(odorTimes(:,2)-odorTimes(:,1)<.5 | odorTimes(:,2)-odorTimes(:,1)>2,:)=[];
    % gather all three regions
    betaData={}; respData={};
    for r=1:length(regions)
        eegData=load(SuperRat(i).([regions{r} 'eegFile']));
        betaPow=eegData.betacontinuous;
        betaPow(:,4)=zscore(betaPow(:,4));
        respPow=eegData.respcontinuous;
        respPow(:,4)=zscore(respPow(:,4));
        clear eegData; % clear this massive file

        % now gather the amplitudes
        [~,~,~,lfpInds,~,times]=event_spikes(betaPow(:,1),odorTimes(:,2),abs(span(1)),span(2));
        betaData{r}=cellfun(@(a) betaPow(a,3), lfpInds,'UniformOutput',false);
        respData{r}=cellfun(@(a) respPow(a,3), lfpInds,'UniformOutput',false);
        % cast all to same size
        for tr=1:length(betaData{r})
            if length(betaData{r}{tr})<dur
                betaData{r}{tr}(end+1:dur)=nan;
                respData{r}{tr}(end+1:dur)=nan;
            elseif length(betaData{r}{tr})>dur
                 betaData{r}{tr}(dur+1:end)=[];
                 respData{r}{tr}(dur+1:end)=[];
            end
        end 
        
        clear betaPow respPow; % dont need this anymore either
    end
    % now get difference in phase
    figure;
    for p=1:3 % pairs
        betaDiff=cellfun(@(a,b) angdiff(a,b), betaData{pairs(p,1)},betaData{pairs(p,2)},'UniformOutput',false);
        subplot(2,3,p);
        betaIM=cell2mat(betaDiff)';
        imagesc(span(1)+1/1500:1/1500:span(2),1:size(betaIM,1),betaIM);
        subplot(2,3,p+3);
        cmean=nan(1,dur); cstd=nan(1,dur);
        for c=1:dur
            cmean(c)=circ_mean(betaIM(:,c));
            cstd(c)=circ_std(betaIM(:,c));
        end
        plot(span(1)+1/1500:1/1500:span(2),cmean); yyaxis right
        plot(span(1)+1/1500:1/1500:span(2),cstd);
        legend({'circ mean','circ std'});
        title(sprintf('%s to %s',regions{pairs(p,1)},regions{pairs(p,2)}))
    end

    clear betaDiff betaIM cmean;
    sgtitle(sprintf('Sess %d Beta',i));    colormap(hsv);
    figure;
    for p=1:3 % pairs
        respDiff=cellfun(@(a,b) angdiff(a,b), respData{pairs(p,1)},respData{pairs(p,2)},'UniformOutput',false);
        subplot(2,3,p);
        respIM=cell2mat(respDiff)';
        imagesc(span(1)+1/1500:1/1500:span(2),1:size(respIM,1),respIM);
        subplot(2,3,p+3);
        cmean=nan(1,dur); cstd=nan(1,dur);
        for c=1:dur
            cmean(c)=circ_mean(respIM(:,c)); 
            cstd(c)=circ_std(respIM(:,c));
        end
        plot(span(1)+1/1500:1/1500:span(2),cmean); yyaxis right
        plot(span(1)+1/1500:1/1500:span(2),cstd);
        legend({'circ mean','circ std'});
        title(sprintf('%s to %s',regions{pairs(p,1)},regions{pairs(p,2)}))
    end
    clear respDiff respIM cmean;
    sgtitle(sprintf('Sess %d RR',i));    colormap(hsv);
end

%% run sessions in aggregate here
% this will get beta offset as a fx of time from odor onset of offset, see
% below to determine wiether its onset or offset

betaOffsets=cell(length(SuperRat),6);
rrOffsets=cell(length(SuperRat),6);
% i think the gist of this is that the beta power just keeps going up the
% longer the animal is in the odor port...
pairs=[1 2; 1 3; 2 3];
span=[-.7 .7]; % secbefore to secafter
dur=1500*diff(span);
for i=1:length(SuperRat)
    % get odor starts
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    odorTimes(odorTimes(:,2)-odorTimes(:,1)<.5 | odorTimes(:,2)-odorTimes(:,1)>2,:)=[];
    % gather all three regions
    betaData={}; respData={};
    for r=1:length(regions)
        eegData=load(SuperRat(i).([regions{r} 'eegFile']));
        betaPow=eegData.betacontinuous;
        betaPow(:,4)=zscore(betaPow(:,4));
        respPow=eegData.respcontinuous;
        respPow(:,4)=zscore(respPow(:,4));
        clear eegData; % clear this massive file

        % now gather the amplitudes
        %
        %  the below line is where you designate locking to onset or offset
        %
        [~,~,~,lfpInds,~,times]=event_spikes(betaPow(:,1),odorTimes(:,2),abs(span(1)),span(2));
        betaData{r}=cellfun(@(a) betaPow(a,3), lfpInds,'UniformOutput',false);
        respData{r}=cellfun(@(a) respPow(a,3), lfpInds,'UniformOutput',false);
        % cast all to same size
        for tr=1:length(betaData{r})
            if length(betaData{r}{tr})<dur
                betaData{r}{tr}(end+1:dur)=nan;
                respData{r}{tr}(end+1:dur)=nan;
            elseif length(betaData{r}{tr})>dur
                 betaData{r}{tr}(dur+1:end)=[];
                 respData{r}{tr}(dur+1:end)=[];
            end
        end 
        
        clear betaPow respPow; % dont need this anymore either
    end
    % now get difference in phase

    for p=1:3 % pairs
        betaDiff=cellfun(@(a,b) angdiff(a,b), betaData{pairs(p,1)},betaData{pairs(p,2)},'UniformOutput',false);
       
        betaIM=cell2mat(betaDiff)';
        cmean=nan(1,dur); cstd=nan(1,dur);
        for c=1:dur
            cmean(c)=circ_mean(betaIM(:,c));
            cstd(c)=circ_std(betaIM(:,c));
        end
        %
        betaOffsets{i,p}=cmean;
        betaOffsets{i,p+3}=cstd;
        
    end

    clear betaDiff betaIM cmean;
    for p=1:3 % pairs
        respDiff=cellfun(@(a,b) angdiff(a,b), respData{pairs(p,1)},respData{pairs(p,2)},'UniformOutput',false);
       
        respIM=cell2mat(respDiff)';
        cmean=nan(1,dur); cstd=nan(1,dur);
        for c=1:dur
            cmean(c)=circ_mean(respIM(:,c)); 
            cstd(c)=circ_std(respIM(:,c));
        end
        rrOffsets{i,p}=cmean;
        rrOffsets{i,p+3}=cstd;
    end
    clear respDiff respIM cmean;
    sgtitle(sprintf('Sess %d RR',i));    colormap(hsv);
end

% for both rhythms is 
% 1 mean PFC-CA1, 
% 2 mean PFC-OB 
% 3 mean CA1-OB 
% 4 std PFC-CA1, 
% 5 std PFC-OB 
% 6 std CA1-OB 
save('E:\Brandeis datasets\Claire Data\ClaireFigs\BetaRROffsetsOdorOn','rrOffsets','betaOffsets','span','pairs','regions');
fprintf('saved Beta and RR offsets to \n %s','E:\Brandeis datasets\Claire Data\ClaireFigs\BetaRROffsetsOdorOff');
%
%
% plotting the above results
%
%
tbins=span(1)+1/1500:1/1500:span(2);
figure(12); figure(22);
for i=1:3
    figure(12);
    subplot(2,3,i);
    betamean=mean(cell2mat(betaOffsets(:,i)),'omitnan');
    betaspread=SEM(cell2mat(betaOffsets(:,i)),1);
    plot(tbins,betamean);
    hold on; patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6);
    title(sprintf('%s,%s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel('Offset ');
    box off; axis tight

    subplot(2,3,i+3);
    betamean=mean(cell2mat(betaOffsets(:,i+3)),'omitnan');
    betaspread=SEM(cell2mat(betaOffsets(:,i+3)),1);
    plot(tbins,betamean);
    hold on; patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6);
    title(sprintf('%s,%s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel('Offset variance');
    box off; axis tight

    figure(22);
    subplot(2,3,i);
    betamean=mean(cell2mat(rrOffsets(:,i)),'omitnan');
    betaspread=SEM(cell2mat(rrOffsets(:,i)),1);
    plot(tbins,betamean);
    hold on; patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6);
    title(sprintf('%s,%s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel('Offset');
    betamean=mean(cell2mat(rrOffsets(:,i+3)),'omitnan');
    box off; axis tight

    subplot(2,3,i+3);
    betaspread=SEM(cell2mat(rrOffsets(:,i+3)),1);
    plot(tbins,betamean);
    hold on; patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6);
    title(sprintf('%s,%s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel('Offset Variance');
    box off; axis tight

end
figure(12); sgtitle('beta'); figure(22); sgtitle('rr');

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
        eegData=load(SuperRat(i).([regions{j} 'eegFile']));
        betaData=eegData.betacontinuous;
        respData=eegData.respcontinuous;
        betaData(:,4)=zscore(betaData(:,4));
        respData(:,4)=zscore(respData(:,4));
        
        mycells=SuperRat(i).units(cellfun(@(a) contains(a,regions{j}), {SuperRat(i).units.area}) &...
            cellfun(@(a) contains(a,celltype), {SuperRat(i).units.type}) & ...
            cellfun(@(a) a(3)<.05, {SuperRat(i).units.taskResponsive}));
        
        % concatenate ALL cells
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
        
        % for session 8, longer trials=later significant spike-beta coherence in ca1
        figure; scatter(odorTimes(:,2)-odorTimes(:,1),phaseTimes);
    end
end
%}       

%% does the beta offset differ between correct trials and incorrect trials?

% this will get beta offset as a fx of time from odor onset of offset, see
% below to determine wiether its onset or offset

betaOffsetsCI=cell(length(SuperRat),6);
rrOffsetsCI=cell(length(SuperRat),6);
% i think the gist of this is that the beta power just keeps going up the
% longer the animal is in the odor port...
pairs=[1 2; 1 3; 2 3];
verbose=0;
for i=1:length(SuperRat)
    % get odor starts
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    % based on data, try last 300 msec before odor release
    
    trialDitch=odorTimes(:,2)-odorTimes(:,1)<.5 | odorTimes(:,2)-odorTimes(:,1)>2;
    odorTimes(:,1)=odorTimes(:,2)-.3;
    % kill trials outside of the proper range
    odorTimes(trialDitch,:)=[];
    corrIncorr10=SuperRat(i).trialdata.CorrIncorr10;
    corrIncorr10(trialDitch,:)=[];
    % gather all three regions
    betaData={}; respData={};
    for r=1:length(regions)
        eegData=load(SuperRat(i).([regions{r} 'eegFile']));
        betaPow=eegData.betacontinuous;
        betaPow(:,4)=zscore(betaPow(:,4)); % unused for now
        respPow=eegData.respcontinuous;
        respPow(:,4)=zscore(respPow(:,4)); % unused for now
        clear eegData; % clear this massive file

        % now gather the amplitudes
        %
        %  the below line is where you designate locking to onset or offset
        %
        [~,~,~,lfpInds,~,times]=event_spikes(betaPow(:,1),odorTimes(:,1),0,diff(odorTimes,1,2));
        betaData{r}=cellfun(@(a) betaPow(a,3), lfpInds,'UniformOutput',false);
        respData{r}=cellfun(@(a) respPow(a,3), lfpInds,'UniformOutput',false);
        % cast all to same size
        %if verbose
        %    mindur=min(cellfun(@(a) length(a), betaData{r}));
        %    betaIM=cell2mat(cellfun(@(a) a(1:mindur), betaData{r},'UniformOutput',false));
        %end
        
        clear betaPow respPow; % dont need this anymore either
    end
    % now get difference in phase

    for p=1:3 % pairs
        betaDiff=cellfun(@(a,b) angdiff(a,b), betaData{pairs(p,1)},betaData{pairs(p,2)},'UniformOutput',false);

        % get mean per trial
        trMean=cell2mat(cellfun(@(a) circ_mean(a), betaDiff,'UniformOutput',false));
        cmean=circ_mean(trMean(corrIncorr10==1));
        imean=circ_mean(trMean(corrIncorr10==0));
        if verbose==1
            mindur=min(cellfun(@(a) length(a), betaDiff));
            betaIM=cell2mat(cellfun(@(a) a(1:mindur), betaDiff,'UniformOutput',false))';
            figure(23);
            subplot(3,3,p);
            imagesc(0:.0015:.3, sum(corrIncorr10==1),betaIM(corrIncorr10==1,:)); 
            title(sprintf('beta C %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
            subplot(3,3,p+3);
            imagesc(0:.0015:.3, sum(corrIncorr10==0),betaIM(corrIncorr10==0,:)); 
            title(sprintf('beta I %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
            subplot(3,3,p+6);
            errorbar([0 1],[imean,cmean],[SEM(trMean(corrIncorr10==0)),SEM(trMean(corrIncorr10==1))])
            hold on; scatter(corrIncorr10,trMean,5,'filled','xjitter','density');
            set(gca,'XTick',[0 1],'XTickLabel',{'Incorrect','correct'})
            colormap hsv
        end
        % incorrect first
        betaOffsetsCI{i,p}=[imean, cmean];
        % maybe need to do a circular dprime here?
        betaOffsetsCI{i,p+3}=dprime(trMean(corrIncorr10==0),trMean(corrIncorr10==1));
        
    end

    clear betaDiff betaIM cmean;
   
    for p=1:3 % pairs
        respDiff=cellfun(@(a,b) angdiff(a,b), respData{pairs(p,1)},respData{pairs(p,2)},'UniformOutput',false);
       
        % get a trial averaged
        trMean=cellfun(@(a) circ_mean(a), respDiff);
        cmean=circ_mean(trMean(corrIncorr10==1));
        imean=circ_mean(trMean(corrIncorr10==0));
        if verbose==1
            mindur=min(cellfun(@(a) length(a), respDiff));
            respIM=cell2mat(cellfun(@(a) a(1:mindur), respDiff,'UniformOutput',false))';
            figure(23);
            subplot(3,3,p);
            imagesc(0:.0015:.3, sum(corrIncorr10==1),respIM(corrIncorr10==1,:)); 
            title(sprintf('rr C %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
            subplot(3,3,p+3);
            imagesc(0:.0015:.3, sum(corrIncorr10==0),respIM(corrIncorr10==0,:)); 
            title(sprintf('rr I %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
            subplot(3,3,p+6);
            errorbar([0 1],[imean,cmean],[SEM(trMean(corrIncorr10==0)),SEM(trMean(corrIncorr10==1))])
            hold on; scatter(corrIncorr10,trMean,5,'filled','xjitter','density');
            set(gca,'XTick',[0 1],'XTickLabel',{'Incorrect','correct'})
            colormap hsv
        end

        % imean first!!!!!!!!
        rrOffsetsCI{i,p}=[imean, cmean];
        rrOffsetsCI{i,p+3}=dprime(trMean(corrIncorr10==0),trMean(corrIncorr10==1));
    end
    clear respDiff respIM cmean;
    sgtitle(sprintf('Sess %d RR',i));    colormap(hsv);
end

% for both rhythms is offset and dprime c-i
% 1 mean PFC-CA1, 
% 2 mean PFC-OB 
% 3 mean CA1-OB 
% 4 dprime PFC-CA1, 
% 5 dprime PFC-OB 
% 6 dprime CA1-OB 
