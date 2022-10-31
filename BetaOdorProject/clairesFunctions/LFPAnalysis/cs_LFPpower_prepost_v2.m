%LFP power pre/post odor

%use filtered LFPs, not PSD

clear
close all
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};

[topDir,figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
freq = 'beta';
for r = 1:length(regions)
    region = regions{r};
    
    allsess = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = [topDir,animal,'Expt\', animal, '_direct\'];
        daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
        tetinfo = loaddatastruct(animDir, animal,'tetinfo');
        odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
        npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
        rewards = loaddatastruct(animDir, animal,'rewards');
        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = daymatrix(daymatrix(:,1) == day,2);
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                %get eeg data
                tet = cs_getMostCellsTet(animal, day, ep, region);
                eeg = loadeegstruct(animDir, animal, 'eeg',day,ep,tet);
                time = geteegtimes(eeg{day}{ep}{tet});
                data = double(eeg{day}{ep}{tet}.data);
                
%                 respfilt = designeegfilt(1500,7,7.5);
%                 filtdata = filtfilt(respfilt, 1, data);
%                 hdata = hilbert(filtdata);
%                 data = abs(hdata);
%                time = eeg{day}{ep}{tet}.starttime + ((0:(length(data(:,1))-1))*(1/eeg{day}{ep}{tet}.samprate));
                
%                 eeg = loadeegstruct(animDir, animal, 'resp',day,ep,tet);
%                 
%                 data = double(eeg{day}{ep}{tet}.data(:,1));
                %data = filtereeg2(
                %get epoch mean for normalization?
%                 switch freq
%                     case 'beta'
                        epmean = mean(data);
                        epsd = std(data);
                    %case 'resp'
%                          rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
%                          rtimes = sort(rtimes);
%                         rewarddata = data(isExcluded(time,rtimes));
%                         epmean = mean(rewarddata);
%                         epsd = std(rewarddata);
                %end
                
                %get pre and post time windows
                nptimes = npWins{day}{ep};
                [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                corrinds = sort([cl;cr]);
                nptimes = nptimes(corrinds,:);
                dur = nptimes(:,2)-nptimes(:,1);
                
                
                
%                 if size(rtimes,1) ~= length(dur)
%                     [rtimes,nptimes,dur] =cs_alignNPReward(rtimes(:,1),nptimes(:,1),dur);
% %                     figure, plot(rtimes(:,1),ones(1,size(rtimes,1)),'r.')
% %                         hold on
% %                         plot(nptimes(:,1),ones(1,size(nptimes,1)),'b.')
%                         %keyboard
%                 end
%                 nptimes = [nptimes, nptimes+dur];
%                
%                 baselinetimes = [rtimes, rtimes + dur];
                
                 baselinetimes = [nptimes(:,1)-dur, nptimes(:,1)];
                
                %get data in windows
                predata = data(isExcluded(time,baselinetimes));
                postdata = data(isExcluded(time,nptimes));
                
%                 %normalize- zscore both to the epoch mean?
%                 
%                 predata = (predata-epmean)/epsd;
%                 postdata = (postdata-epmean)/epsd;
                %normalize post power as a function of pre power?
                
                allsess = [allsess; mean(predata), mean(postdata)];
            end
        end
    end
    load([dataDir, 'PSDdata_',region]);
    post = PSDdata.alldata;
    fpass = PSDdata.fpass;
    stepsize = (fpass(2)-fpass(1))/size(post,2);
    freq = fpass(1):stepsize:fpass(2)-stepsize;
    
    load([dataDir, 'PSDdata_pre_',region]);
    pre = PSDdata.alldata;
    
    mn_pre = mean(pre,1);
    mn_post = mean(post,1);
    figure, plot(freq,mn_pre);
    hold on
    plot(freq,mn_post);
    
    
    CI_pre = [mn_pre-prctile(pre,5); mn_pre + prctile(pre,95)];
    patch([freq,fliplr(freq)],[CI_pre(1,:),fliplr(CI_pre(2,:))],'k','FaceAlpha',0.5)
    
    CI_post = [mn_post-prctile(post,5); mn_pre + prctile(post,95)];
    patch([freq,fliplr(freq)],[CI_post(1,:),fliplr(CI_post(2,:))],'b','FaceAlpha',0.5)
    %set(gca, 'YScale', 'log')
    %
    %      %% RR
    %     RRfreq = find(freq <= 8 & freq>=7);
    %     RR_pre = mean(pre(:,RRfreq),2);
    %     RR_post = mean(post(:,RRfreq),2);
    %
    %     %remove random samples to make sizes equal
    %
    %     diff = length(RR_post)-length(RR_pre);
    %     rem = randsample(length(RR_post),diff);
    %     RR_post(rem) = [];
    %
    %     p = signrank(RR_pre,RR_post);
    %
    %     RRmn_pre = mean(RR_pre);
    %     RRerr_pre = stderr(RR_pre);
    %     RRmn_post = mean(RR_post);
    %     RRerr_post = stderr(RR_post);
    %
    %     figure
    %     hold on
    %     errorbar([1,2],[RRmn_pre, RRmn_post],[RRerr_pre,RRerr_post])
    %      text(1,100,['p = ',num2str(p)]);
    
    
    
    %% Beta
    betafreq = find(freq <= 30 & freq>=15);
    beta_pre = mean(pre(:,betafreq),2);
    beta_post = mean(post(:,betafreq),2);
    
    %remove random samples to make sizes equal
    
    diff = length(beta_post)-length(beta_pre);
    rem = randsample(length(beta_post),diff);
    beta_post(rem) = [];
    
    
    p = signrank(beta_pre,beta_post);
    
    betamn_pre = mean(beta_pre);
    betaerr_pre = stderr(beta_pre);
    betamn_post = mean(beta_post);
    betaerr_post = stderr(beta_post);
    
    errorbar([1,2],[betamn_pre, betamn_post],[betaerr_pre,betaerr_post])
    text(1,100,['p = ',num2str(p)]);
    xticks([1 2])
    xticklabels({'pre','post'});
    xlim([0 3])
    
    figfile = [figDir,'EEG\PowerPrePost_',region];
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
    
    
    
end