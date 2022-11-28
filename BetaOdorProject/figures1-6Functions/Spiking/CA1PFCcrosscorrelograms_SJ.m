% based on predurACCA1GLM4

% SJ: July2015, For Version X6, use to make final figs
clear; %pack;
figdir = '/data25/sjadhav/HPExpt/Figures/RipTrig2015/CrossCorr/';
savefig1=0;
set(0,'defaultaxesfontsize',16);


%% --------------- START: Predicting CA1 SWR activity from preceding AC
%savedirX = '/opt/data15/gideon/HP_ProcessedData/';
savedirX = '/data25/sjadhav/HPExpt/HP_ProcessedData/';

% the following files were created by
% DFSsj_HPexpt_getripalignspiking_ver5.m (I think..)


runscript=0;

if runscript
    
    %load([savedirX 'HP_ripplemod_CA1_alldata_std3_speed4_ntet2_Nspk50_gather_X6'])
    load([savedirX 'HP_ripplemod_CA1_alldata_std3_speed4_ntet2_Nspk50_gather_X6_CA2tag'])
    allripplemodCA1=allripplemod;
    allripplemod_idxCA1=allripplemod_idx;
    
    load([savedirX 'HP_ripplemod_PFC_alldata_std3_speed4_ntet2_Nspk50_gather_X6'])
    %load([savedirX 'HP_ripplemod_PFC_alldata_std3_speed4_ntet2_Nspk50_gather_X6symmetricwindow']);
    
    allripplemodPFC=allripplemod;
    allripplemod_idxPFC=allripplemod_idx;
    
    % combining AC and PFC!!
    % allripplemodAC=[allripplemodAC allripplemodPFC];
    % allripplemod_idxAC=[allripplemod_idxAC ; allripplemod_idxPFC]
    
    combined_idx=unique([allripplemod_idxPFC(:,1:2)],'rows');
    %
    % numtimewindows=size(timewindows,1);
    % ALLSIGRATIOS=NaN(numtimewindows,numtimewindows);
    % ALLERRORRATIOS=NaN(numtimewindows,numtimewindows);
    
    numTrain=1000;
    scrsz = get(0,'ScreenSize');
    % for ACtimewindow=1:numtimewindows
    %     for CA1timewindow=1:numtimewindows
    xx=[-500:500];
    converged1=[];
    n = [];
    nsig = [];
    fracsig=[];
    allErrReal=[];
    allErrShuf=[];
    allPs=[];
    allPsK=[];
    allinds=[];
    allbetas={};
    allsigs={};
    allsigbetasvec=[];
    allcrosscorrs=[];
    allcrosscorrsSWRexc=[];
    allcrosscorrsSWRinh=[];
    allcrosscovsSWRexc=[];
    allcrosscovsSWRinh=[];
    
    plotpairs=0;
    %b=fir1(51,0.05);
    b=gaussian(20,61);
    allpeaks=[];
    allpeaksS=[];
    allpeaksSWRexc=[];
    allpeaksSWRinh=[];
    timer=0;
    
    plotres=0
    
    crosscorrpeaktimes=[];
    pfcsigdifs=[];
    pfcmods=[];
    pfclatepeakinds=[];
    pfcearlypeakinds=[];
    pfclatetroughinds=[];
    pfcearlytroughinds=[];
    pcorrs=[];
    rcorrs=[];
    % swr-exc =1
    % swr-inh =2
    rtimings=[];
    ptimings=[];
    
    % Do Corrlns
    allr=[]; allp=[];
    
    for v=1:size(combined_idx,1)
        curidx=combined_idx(v,[1:2]);%animal day
        PFCidx=find(ismember(allripplemod_idxPFC(:,1:2),curidx,'rows'));
        CA1idx=find(ismember(allripplemod_idxCA1(:,1:2),curidx,'rows'));
        PFCmat=[];
        for j=1:size(PFCidx,1)
            for k=1:size(CA1idx,1)
                PFCind1=PFCidx(j);
                CA1ind1=CA1idx(k);
                % %Orig
                %if allripplemodPFC(1,PFCind1).rasterShufP2<0.05 & allripplemodCA1(1,CA1ind1).rasterShufP2<0.05
                % %X8
                if allripplemodPFC(1,PFCind1).rasterShufP2<0.05 & strcmp(allripplemodPFC(1,PFCind1).FStag,'n') & allripplemodCA1(1,CA1ind1).rasterShufP2<0.05 & strcmp(allripplemodCA1(1,CA1ind1).CA2tag,'n') & strcmp(allripplemodCA1(1,CA1ind1).FStag,'n')
                    PFChist=rast2mat(allripplemodPFC(1,PFCind1).raster);
                    CA1hist=rast2mat(allripplemodCA1(1,CA1ind1).raster);
                    
                    if strcmp(allripplemodCA1(1,CA1ind1).type,'exc')
                        % is the pfc cell swr-excited or swr-inhibited
                        pfcSWRexcinh=strcmp(allripplemodPFC(1,PFCind1).type,'exc');
                        %  mean(mean(PFChist(:,501:650)))>mean(mean(PFChist(:,250:500)));
                        %pfcSWRexcinh=mean(mean(PFChist(:,501:700)))>mean(mean(PFChist(:,50:300)));
                        
                        %now columns are ripples
                        CA1hist2=CA1hist';
                        PFChist2=PFChist';
                        if size(PFChist2,2)==size(CA1hist2,2)
                            PFCh=PFChist2(:);
                            CA1h=CA1hist2(:);
                            currawcrosscorr=xcorr(PFCh,CA1h,500,'coeff')';
                            curcrosscorr=filtfilt(b,1,currawcrosscorr);
                            
                            % this outputs the same as currawcrosscor except the
                            % units, which I think are just co-occurrences
                            spikexcorrOut=spikexcorr(find(PFCh),find(CA1h),1,500);
                            
                            % cross-cov calculation from Siapas paper- verify
                            c1vsc2=spikexcorrOut.c1vsc2;
                            T=length(PFCh)/1000;
                            bin=0.001;
                            spikerate1=sum(PFCh)/T;
                            spikerate2=sum(CA1h)/T;
                            xc = (c1vsc2/(bin*T))-(spikerate1*spikerate2);
                            xcQ= xc*(sqrt((bin*T)./(spikerate1*spikerate2)));
                            nstd=round(0.01/bin); % 10ms/bin
                            g1 = gaussian(nstd, 5*nstd+1);
                            crosscov1 = smoothvect(xcQ, g1);
                            
                            % the PFC cell is SWR-excited
                            if pfcSWRexcinh
                                %allcrosscorrsSWRexc=[allcrosscorrsSWRexc; curcrosscorr];
                                allcrosscovsSWRexc=[allcrosscovsSWRexc; crosscov1];
                                % find peak time
                                [a maxind]=max(curcrosscorr(251:750));
                                allpeaksSWRexc=[allpeaksSWRexc maxind+250];
                                
                            else % the PFC cell is SWR-inhibited
                                %allcrosscorrsSWRinh=[allcrosscorrsSWRinh; curcrosscorr];
                                allcrosscovsSWRinh=[allcrosscovsSWRinh; crosscov1];
                                % find trough time
                                [a minind]=min(curcrosscorr(251:750));
                                allpeaksSWRinh=[allpeaksSWRinh minind+250];
                            end
                            
                            
                            
                            timer=timer+1
                            g2=gaussian(10,51);
                            g3=gaussian(3,3);
                            xaxis=-499:500;
                            % CA1 spiking vs not spiking
                            [firstSpike swrNum]=find(cumsum(cumsum(CA1hist(:,500:700)'))==1);
                            CA1wspikes=CA1hist(swrNum,:);
                            CA1wNospikes=CA1hist(setdiff(1:size(CA1hist,1),swrNum),:);
                            PFCwCA1spikes=PFChist(swrNum,:);
                            PFCwNoCA1spikes=PFChist(setdiff(1:size(CA1hist,1),swrNum),:);
                            pfcampwspikes=sum(PFCwCA1spikes(:,500:700),2);
                            pfcampnospikes=sum(PFCwNoCA1spikes(:,500:700),2);
                            %does the PFC cell fire more on SWRs that the CA1
                            %fired than when the CA1 was silent?
                            pfcsigdif=ranksum(pfcampwspikes,pfcampnospikes);
                            pfcsigdifs=[pfcsigdifs pfcsigdif];
                            % depth of modulation of PFC cell
                            pfcmod=(mean(pfcampwspikes)-mean(pfcampnospikes))/(mean(pfcampwspikes)+mean(pfcampnospikes));
                            pfcmods=[pfcmods pfcmod];
                            
                            % CA1-PFC CrossCorrln
                            CA1nspk = sum(CA1hist(:,500:700),2);
                            PFCnspk = sum(PFChist(:,500:700),2);
                            [rtmp, ptmp] = corrcoef(CA1nspk, PFCnspk);
                            r = rtmp(1,2); p = ptmp(1,2);
                            allr = [allr; r];
                            allp = [allp; p];
                            % CA1 spiking temporal order
                            [sortedFirstSpike sortedFirstSpikeX]=sort(firstSpike);
                            CA1wspikesSorted=CA1wspikes(sortedFirstSpikeX,:);
                            PFCwCA1spikesSorted=PFCwCA1spikes(sortedFirstSpikeX,:);
                            numswrsWspikes=size(CA1wspikes,1);
                            pfcearly=smoothvect(mean(PFCwCA1spikesSorted(1:round(numswrsWspikes/2),:)),g2);
                            pfclate=smoothvect(mean(PFCwCA1spikesSorted(round(numswrsWspikes/2)+1:end,:)),g2);
                            
                            pfcmatsmooth=gaussian(40,50);
                            smoothedpfc=(filtfilt(pfcmatsmooth,1,PFCwCA1spikes')');
                            smoothedpfc=zscore(smoothedpfc')';
                            
                            % For PFC SWR-excited cells, looking at timing only if:
                            % 1. The PFC cell is significantly positively modulated
                            % 2. There are enough SWRs where CA1 fired
                            % 3. The PFC cell fires enough spikes during CA1 firing
                            %%if pfcSWRexcinh & pfcsigdif<0.05 & pfcmod>0 & size(CA1wspikes,1)>20 & sum(sum(PFCwCA1spikes))>30  %orig
                            
                            % Remove condition of significant difference for no spks vs yes spks: pfcsigdiff
                            % %pfcmod>0 is condition for postive corrln for excited cell
                            %if pfcSWRexcinh & pfcmod>0 & size(CA1wspikes,1)>20 & sum(sum(PFCwCA1spikes))>30  %n1 pfcmod condn
                            %if pfcSWRexcinh & size(CA1wspikes,1)>20 & sum(sum(PFCwCA1spikes))>30   %n2 - No condns
                            % Use orig condition of positive corrlns for exc for n3
                            %if pfcSWRexcinh & ~isnan(r) & r>0 & size(CA1wspikes,1)>20 & sum(sum(PFCwCA1spikes))>30 %n3 positive corrln exc condn
                                % Reverse of n3 condtion
                                if pfcSWRexcinh & ~isnan(r) & r<0 & size(CA1wspikes,1)>20 & sum(sum(PFCwCA1spikes))>30  %n3rev negative corrln exc condn
                                
                                [pfclatepeakamp pfclatepeakind]=max(pfclate(500:800));
                                [pfcearlypeakamp pfcearlypeakind]=max(pfcearly(500:800));
                                pfclatepeakinds=[pfclatepeakinds pfclatepeakind];
                                pfcearlypeakinds=[pfcearlypeakinds pfcearlypeakind];
                                [maxamppfc maxamptimepfpc]=max(smoothedpfc(:,500:800)');
                                sigpeaks=find(maxamppfc>1);
                                plotres2=1;
                                %                             if length(sigpeaks>20)
                                %                             [rtiming ptiming]=corrcoef(maxamptimepfpc(sigpeaks),firstSpike(sigpeaks));
                                %                             rtimings=[rtimings rtiming(2,1)];
                                %                             ptimings=[ptimings ptiming(2,1)];
                                %                             end
                                
                                % For PFC SWR-inhibited cells, looking at timing only if:
                                % 1. The PFC cell is significantly negatively modulated
                                % 2. There are enough SWRs where CA1 fired
                                %%elseif ~pfcSWRexcinh & pfcsigdif<0.05 & pfcmod<0 & size(CA1wspikes,1)>20  %orig
                                
                                % Remove condition of significant difference for no spks vs yes spks: pfcsigdiff
                                % pfcmod<0 is condition for negative corrln for inhibited cell
                                %elseif ~pfcSWRexcinh & pfcmod<0 & size(CA1wspikes,1)>20   %n1
                                %%elseif ~pfcSWRexcinh & size(CA1wspikes,1)>20   %n2
                                % Use orig condition of positive corrlns for exc for n3
                            %elseif ~pfcSWRexcinh & ~isnan(r) & r<0 & size(CA1wspikes,1)>20 %n3 negative corrln inh condn
                                % Reverse of n3 condtion
                                elseif ~pfcSWRexcinh & ~isnan(r) & r>0 & size(CA1wspikes,1)>20 %n3rev positive corrln inh condn
                                
                                [pfclatetroughamp pfclatetroughind]=min(pfclate(500:800));
                                [pfcearlytroughamp pfcearlytroughind]=min(pfcearly(500:800));
                                pfclatetroughinds=[pfclatetroughinds pfclatetroughind];
                                pfcearlytroughinds=[pfcearlytroughinds pfcearlytroughind];
                                plotres2=1;
                            else
                                plotres2=0;
                            end
                            
                            
                            if plotres & plotres2
                                if 1%pfcsigdif<0.05 & mean(pfcampwspikes)>mean(pfcampnospikes)&size(CA1wNospikes,1)>30&sum(sum(PFCwCA1spikes))>30
                                    
                                    % CA1 spiking vs not spiking
                                    figure('Position',[30 400 scrsz(3)/4 800]);
                                    g2=gaussian(10,51);
                                    subplot(3,2,1);
                                    imagesc(1-ceil(filtfilt(g3,1,CA1wspikes')'));
                                    title('CA1 w spikes')
                                    subplot(3,2,3);
                                    imagesc(1-ceil(filtfilt(g3,1,CA1wNospikes')'));
                                    colormap(gray)
                                    title('CA1 w No spikes')
                                    subplot(3,2,5);
                                    plot(smoothvect(mean(CA1wspikes),g2),'r');hold on;plot(smoothvect(mean(CA1wNospikes),g2))
                                    legend('CA1 spikes','CA1 No spikes')
                                    subplot(3,2,2);
                                    imagesc(1-ceil(filtfilt(g3,1,PFCwCA1spikes')'));
                                    title('PFC w CA1 spikes')
                                    subplot(3,2,4);
                                    imagesc(1-ceil(filtfilt(g3,1,PFCwNoCA1spikes')'));
                                    colormap(gray)
                                    title('PFC w No CA1 spikes')
                                    subplot(3,2,6)
                                    plot(smoothvect(mean(PFCwCA1spikes),g2),'r');hold on;plot(smoothvect(mean(PFCwNoCA1spikes),g2))
                                    title(pfcsigdif)
                                    
                                    set(0,'defaultAxesFontName', 'Arial')
                                    set(0,'defaultUIControlFontName', 'Arial')
                                    set(0,'defaultTextFontName', 'Arial')
                                    set(0,'defaultaxesfontsize',40);
                                    
                                    figure('Position',[1200 400 scrsz(3)/4 800])
                                    subplot(2,2,1);
                                    imagesc(xaxis,1:size(CA1wspikesSorted,1),1-ceil(filtfilt(g3,1,CA1wspikesSorted')'));colormap(gray);
                                    title('CA1 w spikes sorted')
                                    xlabel('Time (ms)')
                                    ylabel('SWR','FontSize',40)
                                    subplot(2,2,2);
                                    imagesc(xaxis,1:size(PFCwCA1spikesSorted,1),1-ceil(filtfilt(g3,1,PFCwCA1spikesSorted')'));colormap(gray)
                                    title('PFC w CA1 spikes sorted')
                                    xlabel('Time (ms)')
                                    ylabel('SWR')
                                    subplot(2,2,3)
                                    plot(xaxis,1000*smoothvect(mean(CA1wspikesSorted(1:round(numswrsWspikes/2),:)),g2),'linewidth',2)
                                    hold on
                                    plot(xaxis,1000*smoothvect(mean(CA1wspikesSorted(round(numswrsWspikes/2)+1:end,:)),g2),'r','linewidth',2)
                                    legend('early CA1 spikes','late CA1 spikes')
                                    xlabel('Time (ms)')
                                    ylabel('Smoothed mean peri-SWR time histogram (spikes/S)')
                                    subplot(2,2,4)
                                    plot(xaxis,1000*smoothvect(mean(PFCwCA1spikesSorted(1:round(numswrsWspikes/2),:)),g2),'linewidth',2)
                                    hold on
                                    plot(xaxis,1000*smoothvect(mean(PFCwCA1spikesSorted(round(numswrsWspikes/2)+1:end,:)),g2),'r','linewidth',2)
                                    xlabel('Time (ms)')
                                    legend('early CA1 spikes','late CA1 spikes')
                                    if pfcSWRexcinh
                                        title(['erly pk= ' num2str(pfcearlypeakind) ' late pk= ' num2str(pfclatepeakind) ])
                                    else
                                        title(['erly trough= ' num2str(pfcearlytroughind) ' late= ' num2str(pfclatetroughind) ])
                                        
                                    end
                                    ylabel('Smoothed mean peri-SWR time histogram (spikes/S)')
                                    
                                    keyboard
                                    close all
                                end
                            end
                            
                            
                        end
                    end
                end
            end
        end
    end
    
    
    % Save everything for just plotting
    %save([savedirX 'HP_CA1PFC_SWRcrosscorr_X6']); %Orig
    %save([savedirX 'HP_CA1PFC_SWRcrosscorr_X6n1']); %sigdiff condition removed
    %save([savedirX 'HP_CA1PFC_SWRcrosscorr_X6n2']); %sigdiff and pfcmod condition removed
    %save([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n1']); %sigdiff condition removed, no CA2, noFS
    %save([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n2']); %sigdiff and pfcmod condition removed, no CA2, noFS
    %save([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n3_again']); % orig pos corrln condin for exc, neg for inh, again: does it again for confirmation
    save([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n3rev_again']); % reverse of above
else % if runscript
    
    %load([savedirX 'HP_CA1PFC_SWRcrosscorr_X6']); %Orig
    %load([savedirX 'HP_CA1PFC_SWRcrosscorr_X6n1']); %sigdiff condition removed
    %load([savedirX 'HP_CA1PFC_SWRcrosscorr_X6n2']); %sigdiff and pfcmod condition removed
    %load([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n1']); %sigdiff condition removed, no CA2, noFS
    %load([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n2']); %sigdiff and pfcmod condition removed, no CA2, noFS
    %load([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n3_again']); % orig pos corrln condin for exc, neg for inh
    load([savedirX 'HP_CA1PFC_SWRcrosscorr_X8n3rev_again']); % reverse of above
    
end

%%
% xaxis=-500:500;
%
% %figure;imagesc(allcrosscorrs);%caxis([-3 3])
% figure;
% subplot(2,1,1)
% plot(xaxis,filtfilt(b,1,mean(allcrosscorrs)),'k','linewidth',2)
% xlabel('Time (ms)')
% ylabel('CA1-PFC Population-mean non-normalized cross-correlgram')
% subplot(2,1,2)
% plot(xaxis,filtfilt(b,1,mean(allcrosscorrs)),'k','linewidth',2)
% xlabel('Time (ms)')
% ylabel('CA1-PFC Population-mean non-normalized cross-correlgram')
% axis([-150 150 2.05 2.3]);grid;
% %%
% % c=var(allcrosscorrs');
% % [aa bb]=sort(c);
% % figure;imagesc(allcrosscorrs(bb,:)./repmat(mean(allcrosscorrs(bb,1:200),2),1,1001));caxis([0 3])
% xaxis=-499:500;
%
% [cc dd]=sort((mean(allcrosscorrs(:,300:700),2)./mean(allcrosscorrs(:,1:300),2)))
% figure;
%
% imagesc(xaxis,size(allcrosscorrs,1),allcrosscorrs(dd,:)./repmat(mean(allcrosscorrs(dd,1:300),2),1,1001));caxis([0 5])
% xlabel('Time (ms)')
% %% peak timing distributions
% [qreal qq]=hist(allpeaks,250:5:750);
% figure;plot(qq-500,smooth(qreal/sum(qreal)))
% hold on
% [qshuf qq]=hist(allpeaksS,250:5:750);
% plot(qq-500,smooth(qshuf/sum(qshuf)),'r')
% xlabel('Peak time relative to SWR onset (ms)')
% ylabel('fraction of CA1-PFC pairs')


%% cross cov
xaxis=-499:500;
figure;shadedErrorBar(xaxis,mean(allcrosscovsSWRexc,1),std(allcrosscovsSWRexc)./sqrt(size(allcrosscovsSWRexc,1)),'-k',1);
axis([-300 300 -0.2 0.2])
xlabel('Time (ms)')
ylabel('Population standardized cross covariances')
hold on
shadedErrorBar(xaxis,mean(allcrosscovsSWRinh,1),std(allcrosscovsSWRinh)./sqrt(size(allcrosscovsSWRinh,1)),'-g',1);
grid


% Don't shade for making final figure in illustrator
% -------------
figure; hold on;
plot(xaxis,mean(allcrosscovsSWRexc,1),'r-','Linewidth',2);
jbfill(xaxis,mean(allcrosscovsSWRexc,1)+sem(allcrosscovsSWRexc,1),mean(allcrosscovsSWRexc,1)-sem(allcrosscovsSWRexc,1),'r','r',1,1);

plot(xaxis,mean(allcrosscovsSWRinh,1),'b-','Linewidth',2);
jbfill(xaxis,mean(allcrosscovsSWRinh,1)+sem(allcrosscovsSWRinh,1),mean(allcrosscovsSWRinh,1)-sem(allcrosscovsSWRinh,1),'b','b',1,1);

axis([-300 300 -0.2 0.2])
xlabel('Time (ms)')
ylabel('Population standardized cross covariances')
hold on


figfile = [figdir,'RipTrig_CrossCorr_ExcInh']
if savefig1==1,
    print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
end


%% cross cov mean only
xaxis=-499:500;
figure;
plot(xaxis,mean(allcrosscovsSWRexc,1),'color','k','linewidth',2);
axis([-300 300 -0.15 0.15])
xlabel('Time (ms)')
ylabel('Population standardized cross covariances')
hold on
plot(xaxis,mean(allcrosscovsSWRinh,1),'color','g','linewidth',2);
grid
figure;
plot(xaxis,mean(allcrosscovsSWRexc,1),'color','k','linewidth',2);
axis([-80 80 0.06 0.11])
ylim([0.06 0.14]); xlim ([-50 80]);
xlabel('Time (ms)')
ylabel('Population standardized cross covariances')

grid

mean(allpeaksSWRexc), median(allpeaksSWRexc), sem(allpeaksSWRexc)
mean(allpeaksSWRinh), median(allpeaksSWRinh), sem(allpeaksSWRinh)



%%
figure;
subplot(2,1,1)
h=barwitherr([std(pfcearlypeakinds)/sqrt(length(pfcearlypeakinds)) std(pfclatepeakinds)/sqrt(length(pfcearlypeakinds))],[mean(pfcearlypeakinds) mean(pfclatepeakinds)])
set(gca, 'xticklabel',{'early CA1 spiking','late CA1 spiking'})
ylabel('Time of PFC swr-triggered psth peak from swr onset (ms)')
set(h(1),'FaceColor',[0.7 0.7 0.7]);

pexcr=ranksum(pfcearlypeakinds,pfclatepeakinds)
pexc=signrank(pfcearlypeakinds,pfclatepeakinds)

title(['Timing of firing peak of SWR-excited PFC cells, P = ' num2str(pexc)])

subplot(2,1,2)
h=barwitherr([std(pfcearlytroughinds)/sqrt(length(pfcearlytroughinds)) std(pfclatetroughinds)/sqrt(length(pfcearlytroughinds))],[mean(pfcearlytroughinds) mean(pfclatetroughinds)]);
set(gca, 'xticklabel',{'early CA1 spiking','late CA1 spiking'})
ylabel('Time of PFC swr-triggered psth trough from swr onset (ms)')
set(h(1),'FaceColor',[0.7 0.7 0.7]);

pinhr=ranksum(pfcearlytroughinds,pfclatetroughinds)
pinh=signrank(pfcearlytroughinds,pfclatetroughinds)
title(['Timing of firing trough of SWR-inhibited PFC cells, P = ' num2str(pinh)])

[mean(pfcearlypeakinds) mean(pfclatepeakinds)], [sem(pfcearlypeakinds) sem(pfclatepeakinds)]
[mean(pfcearlytroughinds) mean(pfclatetroughinds)], [sem(pfcearlytroughinds) sem(pfclatetroughinds)]


% Make plots for paper fig

figure; hold on;
bar(1, nanmean(pfcearlypeakinds),'r');
bar(2, nanmean(pfclatepeakinds),'r');
errorbar2(1, nanmean(pfcearlypeakinds), nanstderr(pfcearlypeakinds), 0.3, 'r');
errorbar2(2, nanmean(pfclatepeakinds), nanstderr(pfclatepeakinds), 0.3, 'r');
title(['Timing of firing peak of SWR-excited PFC cells, P = ' num2str(pexc)])
ylabel('Time of peak','FontSize',24,'FontWeight','normal')

figdir = '/data25/sjadhav/HPExpt/Figures/RipTrig2015/CrossCorr/';
figfile = [figdir,'CA1PFCExc_earlylatespiking_Rev']
if savefig1==1,
    print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
end


figure; hold on;
bar(1, nanmean(pfcearlytroughinds),'b');
bar(2, nanmean(pfclatetroughinds),'b');
errorbar2(1, nanmean(pfcearlytroughinds), nanstderr(pfcearlytroughinds), 0.3, 'b');
errorbar2(2, nanmean(pfclatetroughinds), nanstderr(pfclatetroughinds), 0.3, 'b');
title(['Timing of firing trough of SWR-inhibited PFC cells, P = ' num2str(pinh)])
ylabel('Time of trough','FontSize',24,'FontWeight','normal')

figdir = '/data25/sjadhav/HPExpt/Figures/RipTrig2015/CrossCorr/';
figfile = [figdir,'CA1PFCInh_earlylatespiking_Rev']
if savefig1==1,
    print('-dpdf', figfile); print('-dpng', figfile, '-r300'); saveas(gcf,figfile,'fig'); print('-depsc2', figfile); print('-djpeg', figfile);
end



figure; hold on; hist(allpeaksSWRexc)
figure; hold on; hist(allpeaksSWRinh)

i=1;




