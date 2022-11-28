%LFP power pre/post odor

%get average PSDs, and also show signed-rank stats for beta and RR

clear
close all
regions = {'CA1','PFC','OB'};

[topDir,figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];

for r = 1:length(regions)
    region = regions{r};
    
    load([dataDir, 'PSDdata_',region]);
    post = PSDdata.sessions_post;
    pre = PSDdata.sessions_pre;
    fpass = PSDdata.fpass;
    stepsize = (fpass(2)-fpass(1))/size(post,2);
    freq = fpass(1):stepsize:fpass(2)-stepsize;
    
    
%     
%     
%     CI_pre = [mn_pre-prctile(pre,5); mn_pre + prctile(pre,95)];
%     patch([freq,fliplr(freq)],[CI_pre(1,:),fliplr(CI_pre(2,:))],'k','FaceAlpha',0.5)
%     
%     CI_post = [mn_post-prctile(post,5); mn_pre + prctile(post,95)];
%     patch([freq,fliplr(freq)],[CI_post(1,:),fliplr(CI_post(2,:))],'b','FaceAlpha',0.5)
     %set(gca, 'YScale', 'log')
%     
%      %% RR 

%use freq bands to exclude theta (8-15)
    RRfreq = find(freq <= 7.95 & freq>=7);
    RR_pre = mean(pre(:,RRfreq),2);
    RR_post = mean(post(:,RRfreq),2);
    
    %remove random samples to make sizes equal
    
    diff = length(RR_post)-length(RR_pre);
    rem = randsample(length(RR_post),diff);
    RR_post(rem) = [];
    
    p = signrank(RR_pre,RR_post);
    
    RRmn_pre = mean(RR_pre);
    RRerr_pre = stderr(RR_pre);
    RRmn_post = mean(RR_post);
    RRerr_post = stderr(RR_post);
   
    figure
    hold on
    errorbar([1,2],[RRmn_pre, RRmn_post],[RRerr_pre,RRerr_post])
     text(1,RRmn_pre,['p = ',num2str(p)]);
     xticks([1 2])
    xticklabels({'pre','post'});
    xlim([0 3])
    
     figfile = [figDir,'EEG\PowerPrePost_RR_',region];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
    
    close
     %% Beta
     %use freq bands to exclude theta (8-15)
    betafreq = find(freq <= 30 & freq>18);
    beta_pre = mean(pre(:,betafreq),2);
    beta_post = mean(post(:,betafreq),2);
    
    
    p = signrank(beta_pre,beta_post);
    
    betamn_pre = mean(beta_pre);
    betaerr_pre = stderr(beta_pre);
    betamn_post = mean(beta_post);
    betaerr_post = stderr(beta_post);
   
    figure
    errorbar([1,2],[betamn_pre, betamn_post],[betaerr_pre,betaerr_post])
    text(1,betamn_pre,['p = ',num2str(p)]);
     xticks([1 2])
    xticklabels({'pre','post'});
    xlim([0 3])
    
figfile = [figDir,'EEG\PowerPrePost_beta_',region];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);

    
close
end