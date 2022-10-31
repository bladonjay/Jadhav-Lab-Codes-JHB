%cs_plotGLMPrediction
clear
close all
[topDir, figDir] = cs_setPaths();

regions = {'CA1','PFC'};
predict = 'Choice';

celltype = 'all';
     switch celltype
        case 'all'
            celltag = '';
        case 'selective'
            celltag = '_selective';
    end

for r = 1:length(regions)
    region = regions{r};
    
    load([topDir, 'AnalysesAcrossAnimals\predict',predict,'GLM_IN_',region,celltag]);
    
    figure
    sem_r = cellfun(@stderr,Acc);
    mn_r = cell2mat(Acc_mn);
    
    errorbar(wins, mn_r, sem_r)
    hold on
    
    %reshape shuffle cell arrays and take mean across sessions. Updated
    %script so that this will not be necessary in the future.
%     test = cellfun(@(x) reshape(x,[100,31]),Acc_shuff,'UniformOutput',false);
%     test2 = cellfun(@(x) mean(x,1),test,'UniformOutput',false);
%     Acc_shuff = cellfun(@(x) x', test2, 'UniformOutput',false);
    
    
    sem_s = cellfun(@stderr,Acc_shuff);
    mn_s = cellfun(@mean,Acc_shuff);
    errorbar(wins, mn_s, sem_s)
    
    labels = [repmat(1,1,length(Acc{1})),repmat(2,1,length(Acc_shuff{1}))];
    data = [cell2mat(Acc);cell2mat(Acc_shuff)];
    
    data_r = cell2mat(Acc);
    %data_r = reshape(data_r,size(data_r,1)*size(data_r,2),1);
    
    data_s = cell2mat(Acc_shuff);
    
    %Multiple rank sum tests with bonferroni 
    alpha = 0.05/size(data_r,2);
    pvals = [];
    for t = 1:size(data_r,2)
        p= ranksum(data_r(:,t),data_s(:,t));
        pvals = [pvals,p];
    end
    
    significant = find(pvals < alpha);
    plot(wins(significant), mn_r(significant)+ sem_r(significant) + 0.02, 'k*');
    
    
    %data_s = reshape(data_s,size(data_s,1)*size(data_s,2),1);
    
%     [p,tbl,stats] = anova2([data_r,data_s],31);
%      multcompare(stats)
%     
%     t = table(labels',data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),...
%         data(:,6),data(:,7),data(:,8),data(:,9),data(:,10),data(:,11),...
%         data(:,12),data(:,13),data(:,14),data(:,15),data(:,16),data(:,17),...
%         data(:,18),data(:,19),data(:,20),'VariableNames',{'DataType','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14','t15','t16','t17','t18','t19','t20'});
%     rm = fitrm(t,'t1-t20~DataType');
%     
%     ranovatbl = ranova(rm)
%     multcompare(ranovatbl)
    
%     figure,
%     %plot(bins,mn,'LineWidth',3);
%     hold on
%     patch([bins,fliplr(bins)], [shuffmn+shuffstd, shuffmn-shuffstd],'k','FaceAlpha',0.3);
%     plot([0 bins(end)],[0.5 0.5],'k--');
    ylabel('Fraction of trials correctly predicted');
    xlabel('Time Window Used for Prediction');
    
    figfile = [figDir,'Interneurons\GLMprediction_IN_',predict,'_',region,celltag];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
end