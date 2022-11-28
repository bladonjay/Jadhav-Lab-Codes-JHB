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


regions={'PFC','CA1','OB'}; pairs=[1 2; 1 3; 2 3];
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
span=[-.7 .7]; % secbefore to secafter
dur=1500*diff(span);
for i=[2 4 9 13 16 23 30 35] %length(SuperRat)
    % get odor starts
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialditch=odorTimes(:,2)-odorTimes(:,1)<.5 | odorTimes(:,2)-odorTimes(:,1)>2;
    odorTimes(trialditch,:)=[];
    corrIncorr10=SuperRat(i).trialdata.CorrIncorr10(trialditch==0);
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
        %     Specify here if you want to align to start or end of sampling
        %
        [~,~,~,windowInds,~,times]=event_spikes(betaPow(:,1),odorTimes(:,2),abs(span(1)),span(2));
        betaData{r}=cellfun(@(a) betaPow(a,3), windowInds,'UniformOutput',false);
        respData{r}=cellfun(@(a) respPow(a,3), windowInds,'UniformOutput',false);
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
        phaseDiff=cellfun(@(a,b) angdiff(a,b), betaData{pairs(p,1)},betaData{pairs(p,2)},'UniformOutput',false);
        subplot(2,3,p);
        phaseMat=cell2mat(phaseDiff)';
        imagesc(span(1)+1/1500:1/1500:span(2),1:size(phaseMat,1),phaseMat);
        subplot(2,3,p+3);
        timeMean=nan(1,dur); cstd=nan(1,dur);
        for c=1:dur
            timeMean(c)=circ_mean(phaseMat(:,c));
            cstd(c)=circ_std(phaseMat(:,c));
        end
        plot(span(1)+1/1500:1/1500:span(2),timeMean); yyaxis right
        plot(span(1)+1/1500:1/1500:span(2),cstd);
        legend({'circ mean','circ std'});
        title(sprintf('%s to %s',regions{pairs(p,1)},regions{pairs(p,2)}))
    end

    clear phaseDiff phaseMat timeMean;
    sgtitle(sprintf('Sess %d Beta',i));    colormap(hsv);
    %{
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
    %}
end

%% run sessions in aggregate here
%
% is beta/ rr phase offset reliably across sessions
%
%
% this will get beta offset as a fx of time from odor onset of offset, see
% below to determine wiether its onset or offset

betaOffsets=cell(length(SuperRat),6);
rrOffsets=cell(length(SuperRat),6);
% i think the gist of this is that the beta power just keeps going up the
% longer the animal is in the odor port...
span=[-1.5 1.5]; % secbefore to secafter
dur=1500*diff(span);
for i=1:length(SuperRat)
    tic;
    % get odor starts
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend, ...
        SuperRat(i).trialdata.CorrIncorr10];
    % use only correct trials
    % kill trials that are too long or too short
    odorTimes(odorTimes(:,2)-odorTimes(:,1)<.5 | odorTimes(:,2)-odorTimes(:,1)>2,:)=[];
    % gather all three regions
    betaData={}; respData={};
    for r=1:length(regions)
        eegData=load(SuperRat(i).([regions{r} 'eegFile']));
        betaPow=eegData.betacontinuous;
        %betaPow(:,4)=zscore(betaPow(:,4));
        %betaPow(betaPow(:,4)<2,3)=nan;
        respPow=eegData.respcontinuous;
        %respPow(:,4)=zscore(respPow(:,4));
        %respPow(respPow(:,4)<2,3)=nan;
        clear eegData; % clear this massive file

        % now gather the amplitudes
        %
        %  the below line is where you designate locking to onset or offset
        %
        [~,~,~,windowInds]=event_spikes(betaPow(:,1),odorTimes(:,2),abs(span(1)),span(2));
        %[~,~,~,odorInds]=event_spikes(betaPow(:,1),odorTimes(:,1),0,odorTimes(:,2)-odorTimes(:,1));
        [~,~,~,odorInds]=event_spikes(betaPow(:,1),odorTimes(:,2),.5,0);
        

        % get the phase of the beta and rr
        betaData{r}=cellfun(@(a) betaPow(a,3), windowInds,'UniformOutput',false);
        respData{r}=cellfun(@(a) respPow(a,3), windowInds,'UniformOutput',false);
        betaData{r+3}=cellfun(@(a) betaPow(a,3), odorInds,'UniformOutput',false);
        respData{r+3}=cellfun(@(a) respPow(a,3), odorInds,'UniformOutput',false);

        % cast all to same size for plotting
        for tr=1:length(betaData{r})
            if length(betaData{r}{tr})<dur
                betaData{r}{tr}(end+1:dur)=nan;
                respData{r}{tr}(end+1:dur)=nan;
            elseif length(betaData{r}{tr})>dur
                betaData{r}{tr}(dur+1:end)=[];
                respData{r}{tr}(dur+1:end)=[];
            end
        end
        trkill=cellfun(@(a) isempty(a), betaData{r+3});
        
        clear betaPow respPow; % dont need this anymore either
    end
    % now get difference in phase
    % betaData is cell array of ntrials, m timebins

    for p=1:3 % pairs
        % for each pair calculate the offset for each trial and each
        % timebin, and also a shuffled answer (say 500 times)
        phaseDiff=cellfun(@(a,b) angdiff(a,b), betaData{pairs(p,1)},betaData{pairs(p,2)},'UniformOutput',false);
        phaseMat=cell2mat(phaseDiff)';

        % mean each trial across time (for trial need to pull odor period)
        % this will squish out empty cells, so you have to keep track below
        trDiff=cell2mat(cellfun(@(a,b) circ_mean(angdiff(a,b)),...
            betaData{pairs(p,1)+3},betaData{pairs(p,2)+3},'UniformOutput',false));
        
        % mean across trials
        timeMean=nan(1,dur); 
        % circ_mean only works on vectors
        for c=1:dur
            timeMean(c)=circ_mean(phaseMat(:,c));
        end

        % tabulate
        betaOffsets{i,p}=timeMean; % 1 x 2100 vector, for the true offsets across trials
        betaOffsets{i,p+3}=circ_mean(trDiff); % all trials, but only odor period
        betaOffsets{i,p+6}=trDiff(odorTimes(:,3)==1 & ~trkill'); % correct trials
        betaOffsets{i,p+9}=trDiff(odorTimes(:,3)==0 & ~trkill'); % incorrect tirlas

        clear timeMean;

        phaseDiff=cellfun(@(a,b) angdiff(a,b), respData{pairs(p,1)},respData{pairs(p,2)},'UniformOutput',false);
        respMAT=cell2mat(phaseDiff)';
        trDiff=cell2mat(cellfun(@(a,b) circ_mean(angdiff(a,b)),...
            respData{pairs(p,1)},respData{pairs(p,2)},'UniformOutput',false));

        timeMean=nan(1,dur);
        for c=1:dur
            timeMean(c)=circ_mean(respMAT(:,c));  % mean is still circular
        end
        rrOffsets{i,p}=timeMean;
        rrOffsets{i,p+3}=circ_mean(trDiff); % all trials, but only odor period
        rrOffsets{i,p+6}=trDiff(odorTimes(:,3)==1 & ~trkill'); % correct trials
        rrOffsets{i,p+9}=trDiff(odorTimes(:,3)==0 & ~trkill'); % incorrect tirlas
    end
    clear phaseDiff respMAT timeMean trDiff;
    fprintf('running session %d, this sess took %.2f mins \n',i,toc/60)
end

% output will be rrOffsets betaOffsets, each is a n session by 9
% for both rhythms is 
% 1 mean PFC-CA1, -1.5 to 1.5
% 2 mean PFC-OB 
% 3 mean CA1-OB 
% 4 PFC-CA1 correct (500msec)
% 5 PFC-OB 
% 6 CA1-OB 
% 7 PFC-CA1 incorrect
save('E:\Brandeis datasets\Claire Data\ClaireFigs\BetaRROffsetsOdorOff5','rrOffsets','betaOffsets','span','pairs','regions');
fprintf('saved Beta and RR offsets to \n %s','E:\Brandeis datasets\Claire Data\ClaireFigs\BetaRROffsetsOdorOff5');
%%
load('E:\Brandeis datasets\Claire Data\ClaireFigs\BetaRROffsetsOdorOff5.mat')
%%
%
% plotting the above results
%
%
tbins=span(1)+1/1500:1/1500:span(2);
figure(12);
for i=1:3
    figure(12);
    subplot(2,3,i);
    betamean=[]; betaspread=[];
    betamat=cell2mat(betaOffsets(:,i));
    for ts=1:size(betamat,2)
        betamean(ts)=circ_mean(betamat(:,ts));
        betaspread(ts)=circ_confmean(betamat(:,ts),.01);
    end
    plot(tbins,cell2mat(betaOffsets(:,i)),'Color',[.8 .8 .8]); hold on;
    plot(tbins,betamean,'Color',colors(mod(i,3)+1,:),'LineWidth',2);
    patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6); 
    title(sprintf('%s,%s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel(sprintf('Phase Offset \n %s <-> %s', regions{pairs(i,1)},regions{pairs(i,2)}));

    box off; axis tight
    plot(get(gca,'XLim'),[0 0],'k')
    % we'll do this trialwise
    
    betasess=[]; betaspread=[];
    for ses=1:size(betaOffsets,1) % for each session
        betasess(ses)=circ_mean(betaOffsets{ses,i}(tbins>-.5 & tbins <0));
    end
    [pval,m,upper,lower]=circ_mtest(betasess,0);
    title(sprintf('beta x=%.2f  (%.2f to %.2f) Ha %d',m,lower,upper,pval))

     subplot(2,3,i+3);

    betamean=[]; betaspread=[];
    betamat=cell2mat(rrOffsets(:,i));
    for ts=1:size(betamat,2)
        betamean(ts)=circ_mean(betamat(:,ts));
        betaspread(ts)=circ_confmean(betamat(:,ts),.01);
    end
    plot(tbins,cell2mat(rrOffsets(:,i)),'Color',[.8 .8 .8]); hold on;

    plot(tbins,betamean,'Color',colors(mod(i,3)+1,:),'LineWidth',2);
    patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6);
    title(sprintf('%s - %s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel(sprintf('Phase Offset \n %s <-> %s', regions{pairs(i,1)},regions{pairs(i,2)}));
    box off; axis tight
    plot(get(gca,'XLim'),[0 0],'k');

    betasess=[]; betaspread=[];
    for ses=1:size(rrOffsets,1) % for each session
        betasess(ses)=circ_mean(rrOffsets{ses,i}(tbins>-.5 & tbins <0));
    end
    [pval,m,upper,lower]=circ_mtest(betasess,0);
    title(sprintf('beta x=%.2f  (%.2f to %.2f) Ha %d',m,lower,upper,pval))

end
figure(12); sgtitle('beta top RR bottom'); 
%% report p values for phase offset as compared to zero
% * THIS IS DURING ODOR SAMPLING EPOCH ONLY
nboots=100000;
figure; plot([0.5 3.5],[0 0],'k'); hold on;
betaRRstats=[];
for i=1:3
    % lets get the mean offset for correct and incorrect trials
    betamat=cell2mat(betaOffsets(:,i+3)); % all trials
    scatter(repmat(i-.1,1,length(betamat)), betamat,8,rhythmcolors(1,:),'filled',...
                'XJitter','density','XJitterWidth',.2,'MarkerFaceAlpha',.4)
    betanull=nan(nboots,1);
    randns=(randi(2,size(betamat,1),nboots)-1.5)*2;
    nullmat=betamat.*randns;
    parfor bt=1:nboots
        betanull(bt)=circ_mean(nullmat(:,bt))
    end
    pval=min([mean(circ_mean(betamat)>betanull) mean(circ_mean(betamat)<betanull)]);
    fprintf('Beta offset %s to %s = %.2f, p=%.2e \n',regions{pairs(i,1)},regions{pairs(i,2)},...
        circ_mean(betamat), pval);

    rrmat=cell2mat(rrOffsets(:,i+3)); % all trials
    scatter(repmat(i+.1,1,length(rrmat)), rrmat,8,rhythmcolors(2,:),'filled',...
        'XJitter','density','XJitterWidth',.2,'MarkerFaceAlpha',.4)
    rrnull=nan(nboots,1);
    nullmat=rrmat.*randns;
    parfor bt=1:nboots
        rrnull(bt)=circ_mean(nullmat(:,bt))
    end
    pval=min([mean(circ_mean(rrmat)>rrnull) mean(circ_mean(rrmat)<rrnull)]);
    fprintf('RR offset %s to %s = %.2f, p=%.2e \n',regions{pairs(i,1)},regions{pairs(i,2)},...
        circ_mean(rrmat), pval);
    betaRRstats(i,1)=circ_mean(betamat);
    betaRRstats(i,2)=circ_std(betamat);
    betaRRstats(i,3)=circ_mean(rrmat);
    betaRRstats(i,4)=circ_std(rrmat);
        clear betanull rrnull betamat nullmat

end
hold on;
eb=errorbar([1:3]-.1, betaRRstats(:,1),betaRRstats(:,2),'.','CapSize',0,'LineWidth',2, ...
    'color',rhythmcolors(1,:),'MarkerSize',15)
eb(2)=errorbar([1:3]+.1, betaRRstats(:,3),betaRRstats(:,4),'.','CapSize',0,'LineWidth',2,...
    'color',rhythmcolors(2,:),'MarkerSize',15)
set(gca,'YLim',[-pi pi],'xlim',[.5 3.5]);
set(gca,'XTick',[1:3],'XTickLabel',{'PFC-CA1','PFC-OB','CA1-OB'});
ylabel('Offset in Radians');
legend(eb,'Beta Rhythm','RR'); legend('boxoff');
box off;
%% now see correct vs incorrect
% run sessionwise dprime, report mean effect size, and then report ##
% passing criterion.  then you can bootstrap mean effect size, and ##
% passing criterion (p=0.05)
% the bootstrap can be a shuffle of c and I trials by catting, randperm i
% and taking first ncorr trials
pretext={'correct','Incorect'};
nboots=10000;

    for i=1:3
        % for each pair, get the null by randomly signing each offset
        % null
        % first correct, then incorrect
        betamat=cellfun(@(a,b) angdiff(circ_mean(a), circ_mean(b)), betaOffsets(:,i+6), betaOffsets(:,i+9)); % all trials
        betanull=nan(nboots,1);
        randns=(randi(2,size(betamat,1),nboots)-1.5)*2;
        nullmat=betamat.*randns;
        parfor bt=1:nboots
            betanull(bt)=circ_mean(nullmat(:,bt))
        end
        pval=min([mean(circ_mean(betamat)>betanull) mean(circ_mean(betamat)<betanull)]);
        fprintf('Beta offset %s to %s = %.2f, p=%.2e \n',regions{pairs(i,1)},regions{pairs(i,2)},...
            circ_mean(betamat), pval);

        rrmat=cell2mat(rrOffsets(:,i+3)); % all trials
        rrnull=nan(nboots,1);
        nullmat=rrmat.*randns;
        parfor bt=1:nboots
            rrnull(bt)=circ_mean(nullmat(:,bt))
        end
        pval=min([mean(circ_mean(rrmat)>rrnull) mean(circ_mean(rrmat)<rrnull)]);
        fprintf('%s RR offset %s to %s = %.2f, p=%.2e \n',pretext{k},regions{pairs(i,1)},regions{pairs(i,2)},...
            circ_mean(rrmat), pval);
        clear betanull rrnull betamat nullmat
    end



%% this is a close look at the pvalue at each timepoint
%
%  grab significance with a bootstrap
%

nboots=1000;
% plotting the results with P value on yyaxis right
tbins=span(1)+1/1500:1/1500:span(2);
figure(12);
for i=1:3
    figure;
    subplot(2,1,1)
    betamean=[]; betaspread=[];
    betamat=cell2mat(betaOffsets(:,i));
    for ts=1:size(betamat,2)
        betamean(ts)=circ_mean(betamat(:,ts));
        betaspread(ts)=circ_confmean(betamat(:,ts),.01);
    end
    plot(tbins,cell2mat(betaOffsets(:,i)),'Color',[.8 .8 .8]); hold on;
    plot(tbins,betamean,'Color',colors(mod(i,3)+1,:),'LineWidth',2);
    patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6); 
    title(sprintf('%s,%s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel(sprintf('Phase Offset \n %s <-> %s', regions{pairs(i,1)},regions{pairs(i,2)}));

    box off; axis tight
    plot(get(gca,'XLim'),[0 0],'k')
    % we'll do this trialwise

    %[pval,m,upper,lower]=circ_mtest(betamean,0);
    %title(sprintf('beta x=%.2f  (%.2f to %.2f) Ha %d',m,lower,upper,pval))

    % null
    for bt=1:nboots
        randns=(randi(2,size(betamat,1),size(betamat,2))-1.5)*2;
        bootmat=betamat.*randns;
        for ts=1:size(betamat,2)
            betanull(bt,ts)=circ_mean(betamat(:,ts).*randns(:,ts));
        end
    end
    yyaxis right;
    plot(tbins,mean(betamean<betanull),'k'); set(gca,'YScale','log');
    ylabel('P value');
    
    
    subplot(2,1,2)

    betamean=[]; betaspread=[];
    betamat=cell2mat(rrOffsets(:,i));
    for ts=1:size(betamat,2)
        betamean(ts)=circ_mean(betamat(:,ts));
        betaspread(ts)=circ_confmean(betamat(:,ts),.01);
    end
    plot(tbins,cell2mat(rrOffsets(:,i)),'Color',[.8 .8 .8]); hold on;

    plot(tbins,betamean,'Color',colors(mod(i,3)+1,:),'LineWidth',2);
    patch([tbins fliplr(tbins)],[betamean+betaspread fliplr(betamean-betaspread)],...
        colors(mod(i,3)+1,:),'LineStyle','none','FaceAlpha',.6);
    title(sprintf('%s - %s',regions{pairs(i,1)},regions{pairs(i,2)}));
    xlabel('time from odor offset'); ylabel(sprintf('Phase Offset \n %s <-> %s', regions{pairs(i,1)},regions{pairs(i,2)}));
    box off; axis tight
    plot(get(gca,'XLim'),[0 0],'k');

    [pval,m,upper,lower]=circ_mtest(betamean,0);
    title(sprintf('beta x=%.2f  (%.2f to %.2f) Ha %d',m,lower,upper,pval))
    % null
    for bt=1:nboots
        randns=(randi(2,size(betamat,1),size(betamat,2))-1.5)*2;
        bootmat=betamat.*randns;
        for ts=1:size(betamat,2)
            betanull(bt,ts)=circ_mean(betamat(:,ts).*randns(:,ts));
        end
    end
    yyaxis right;
    plot(tbins,mean(betamean<betanull)); set(gca,'YScale','log');
    ylabel('P value');


end

%% this is some workspace here
figure(12); sgtitle('beta top RR bottom'); 

% bootstrap now for significance
% for the bootstrap i just need to randomize the sign of each betaoffset
% value 500 times
rstream = RandStream('dsfmt19937','Seed',16);
RandStream.setGlobalStream(rstream);

tbins=span(1)+1/1500:1/1500:span(2);
figure(12);
nboots=1000;
for i=1:3

    betamean=[]; betanull=[];
    % this is the grey lines in the plot above, now we randomize sign of
    % those
    betamat=cell2mat(betaOffsets(:,i));
    % circ_mean doesnt handle matrices
    for bt=1:nboots
        randns=(randi(2,size(betamat,1),size(betamat,2))-1.5)*2;
        bootmat=betamat.*randns;
        for ts=1:size(betamat,2)

            betanull(bt,ts)=circ_mean(betamat(:,ts).*randns(:,ts));
        end
    end
    for ts=1:size(betamat,2)
        betamean(ts)=circ_mean(betamat(:,ts));
    end
    % so my question is for -.5 to 0 is the grand average significant
    figure; plot(mean(betamean<betanull));
end


%% so we know that beta power increases with time, but what about spike-beta coherence
%{

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
%
%
%  need to formalize this and report the results in the paper
%
%
%
% this will get beta offset as a fx of time from odor onset of offset, see
% below to determine wiether its onset or offset
bytrial=0;

if bytrial==1
    betaOffsetsCI=cell(length(SuperRat),6);
    rrOffsetsCI=cell(length(SuperRat),6);
else
    betaOffsetsCI=cell(3,2);
    rrOffsetsCI=cell(3,2);
end

% i think the gist of this is that the beta power just keeps going up the
% longer the animal is in the odor port...
pairs=[1 2; 1 3; 2 3];
verbose=0;
for i=1:length(SuperRat)
    % get odor starts
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    % based on data, try last 300 msec before odor release
    
    trialDitch=odorTimes(:,2)-odorTimes(:,1)<.5 | odorTimes(:,2)-odorTimes(:,1)>2;

    % CS42 has last trial cut off in data for every session...
    %odorTimes(:,1)=odorTimes(:,2)-.5;
    % kill trials outside of the proper range
    odorTimes(trialDitch,:)=[];
    corrIncorr10=SuperRat(i).trialdata.CorrIncorr10;
    corrIncorr10(trialDitch,:)=[];
    % gather all three regions
    betaData={}; respData={}; trialkill=zeros(length(odorTimes),1);
    for r=1:length(regions)
        eegData=load(SuperRat(i).([regions{r} 'eegFile']));
        betaPow=eegData.betacontinuous;     
        respPow=eegData.respcontinuous;
        clear eegData; % clear this massive file

        % now gather the amplitudes
        %
        %  the below line is where you designate locking to onset or offset
        %
        %[~,~,~,lfpInds,~,times]=event_spikes(betaPow(:,1),odorTimes(:,1),0,diff(odorTimes,1,2));

        [~,~,~,windowInds,~,times]=event_spikes(betaPow(:,1),odorTimes(:,2),1,.5);
        newx=-1:(1/1500):.5-(1/1500);

        betaData{r}=cellfun(@(a) betaPow(a,3), windowInds,'UniformOutput',false);
        respData{r}=cellfun(@(a) respPow(a,3), windowInds,'UniformOutput',false);

        % some sessions have missing data, this removes those missing
        % points
        trialkill(cellfun(@(a) isempty(a), betaData{r}))=1;
        trialkill(cellfun(@(a) isempty(a), respData{r}))=1;
        
        clear betaPow respPow; % dont need this anymore either
    end
        
    % now standardize the sizes of these
    for r=1:3
        betaData{r}(trialkill==1)=[];
        respData{r}(trialkill==1)=[];
        rightdur=mode(cellfun(@(a) length(a), betaData{r}));
        betaData{r}=cellfun(@(a) fixMatSize(a,rightdur), betaData{r},'UniformOutput',false);
        respData{r}=cellfun(@(a) fixMatSize(a,rightdur), respData{r},'UniformOutput',false);

    end
    corrIncorr10(trialkill==1)=[];
    % now get difference in phase

    for p=1:3 % pairs
        phaseDiff=cellfun(@(a,b) angdiff(a,b), betaData{pairs(p,1)},betaData{pairs(p,2)},'UniformOutput',false);

        % get mean per trial 
        if bytrial==1
            trMean=cell2mat(cellfun(@(a) circ_mean(a), phaseDiff,'UniformOutput',false));
            timeMean=circ_mean(trMean(corrIncorr10==1));
            imean=circ_mean(trMean(corrIncorr10==0));
            betaOffsetsCI{i,p}=[imean, timeMean];
            betaOffsetsCI{i,p+3}=dprime(trMean(corrIncorr10==0),trMean(corrIncorr10==1));
        else
            %The alternative is to build a time matrix...
            trMat=cell2mat(phaseDiff)';
            timeMean=nan(size(trMean,2),1); imean=timeMean;
            for tm=1:size(trMat,2)
                timeMean(tm)=circ_mean(trMat(corrIncorr10==1,tm));
                imean(tm)=circ_mean(trMat(corrIncorr10==0,tm));
            end
            betaOffsetsCI{p,1}(i,:)=timeMean;
            betaOffsetsCI{p,2}(i,:)=imean;
        end

        if verbose==1
            if bytrial==1
            mindur=min(cellfun(@(a) length(a), phaseDiff));
            phaseMat=cell2mat(cellfun(@(a) a(1:mindur), phaseDiff,'UniformOutput',false))';
            figure(23);
            subplot(3,3,p);
            imagesc(0:.0015:.3, sum(corrIncorr10==1),phaseMat(corrIncorr10==1,:)); 
            title(sprintf('beta C %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
            subplot(3,3,p+3);
            imagesc(0:.0015:.3, sum(corrIncorr10==0),phaseMat(corrIncorr10==0,:)); 
            title(sprintf('beta I %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
            subplot(3,3,p+6);
            errorbar([0 1],[imean,timeMean],[SEM(trMean(corrIncorr10==0)),SEM(trMean(corrIncorr10==1))])
            hold on; scatter(corrIncorr10,trMean,5,'filled','xjitter','density');
            set(gca,'XTick',[0 1],'XTickLabel',{'Incorrect','correct'})
            colormap hsv
            else
                figure; 
                subplot(2,2,1); imagesc(newx,[],trMat(corrIncorr10==1,:))
                subplot(2,2,2); imagesc(newx,[],trMat(corrIncorr10==0,:));
                s=subplot(2,2,3); plot(newx,timeMean(:,1));
                s(2)=subplot(2,2,4); plot(newx,imean(:,1));
                colormap hsv
                linkaxes(s);
            end
        end
        % incorrect first
        
        % maybe need to do a circular dprime here?
        %
        
    end

    clear phaseDiff phaseMat timeMean;

%     for p=1:3 % pairs
%         respDiff=cellfun(@(a,b) angdiff(a,b), respData{pairs(p,1)},respData{pairs(p,2)},'UniformOutput',false);
% 
%         % get a trial averaged
%         trMean=cellfun(@(a) circ_mean(a), respDiff);
%         cmean=circ_mean(trMean(corrIncorr10==1));
%         imean=circ_mean(trMean(corrIncorr10==0));
%         if verbose==1 && bytrial==1
%             mindur=min(cellfun(@(a) length(a), respDiff));
%             respIM=cell2mat(cellfun(@(a) a(1:mindur), respDiff,'UniformOutput',false))';
%             figure(23);
%             subplot(3,3,p);
%             imagesc(0:.0015:.3, sum(corrIncorr10==1),respIM(corrIncorr10==1,:));
%             title(sprintf('rr C %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
%             subplot(3,3,p+3);
%             imagesc(0:.0015:.3, sum(corrIncorr10==0),respIM(corrIncorr10==0,:));
%             title(sprintf('rr I %s vs %s',regions{pairs(p,1)},regions{pairs(p,2)}));
%             subplot(3,3,p+6);
%             errorbar([0 1],[imean,cmean],[SEM(trMean(corrIncorr10==0)),SEM(trMean(corrIncorr10==1))])
%             hold on; scatter(corrIncorr10,trMean,5,'filled','xjitter','density');
%             set(gca,'XTick',[0 1],'XTickLabel',{'Incorrect','correct'})
%             colormap hsv
%         end
% 
%         % imean first!!!!!!!!
%         rrOffsetsCI{i,p}=[imean, cmean];
%         rrOffsetsCI{i,p+3}=dprime(trMean(corrIncorr10==0),trMean(corrIncorr10==1));
%     end
%     clear respDiff respIM cmean;
%     sgtitle(sprintf('Sess %d RR',i));    colormap(hsv);
end

% for both rhythms is offset and dprime c-i
% 1 mean PFC-CA1, 
% 2 mean PFC-OB 
% 3 mean CA1-OB 
% 4 dprime PFC-CA1, 
% 5 dprime PFC-OB 
% 6 dprime CA1-OB 
%%
% now plot this out as a paired plot
if bytrial==1
% for each pair
for p=1:3
    subplot(1,3,p);
    ICoffsets=cell2mat(betaOffsetsCI(:,p));
    plot([1 2],ICoffsets,'Color',[.8 .8 .8]); hold on;

    errorbar([1 2], mean(ICoffsets),SEM(ICoffsets),'LineWidth',2)
    xlim([.9 2.1])
    ranksum(ICoffsets(:,1),ICoffsets(:,2))
    circ_cmtest(ICoffsets(:), linearize(repmat([1 2],length(ICoffsets),1)))
end
end
% it looks like the phase offsets are actually pretty tightly locked to the
% event, so i may have to do C-I by time locked to offset

if bytrial==0
    for p=1:3
        figure;
        subplot(1,2,1)
        imagesc(betaOffsetsCI{p,1});
        subplot(1,2,2);
        imagesc(betaOffsetsCI{p,2});
        figure;
        plot(angdiff(betaOffsetsCI{p,1},betaOffsetsCI{p,2})');
    end
end