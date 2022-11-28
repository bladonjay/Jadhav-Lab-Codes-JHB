[topDir, figDir] = cs_setPaths();
regions = {'CA1','PFC'};

load([topDir,'AnalysesAcrossAnimals\singleCellSelectivity.mat']);
close all


figure, hold on
colors = {'b','r'};
for r = 1:length(regions)
    region = regions{r};
    selectivityTimes = vertcat(singleCellSelectivity.(region){:,2});
    h = histogram(selectivityTimes,[0:0.1:1],'facealpha',0.5);
    plot([mean(selectivityTimes) mean(selectivityTimes)], [0 8], [colors{r},'--'],'LineWidth',2.5)
end




% allTimes = []; grp = [];
% for r = 1:length(regions)
%     region = regions{r};
%     selectivityTimes = vertcat(singleCellSelectivity.(region){:,2});
%     grp= [grp; repmat(r,size(selectivityTimes,1),1)];
%     allTimes = vertcat(allTimes, selectivityTimes);
% end
% figure, hold on
% boxplot(allTimes, grp,'symbol','k.','labels',regions)
% %xlabel(regions)
% ylabel('Time of selectivity emergence (seconds)')
% jitter = -0.03+rand(length(grp),1)*(0.06);
% grp = grp+jitter;
% plot(grp,allTimes,'k.');

[~,p] = ttest2(vertcat(singleCellSelectivity.CA1{:,2}), vertcat(singleCellSelectivity.PFC{:,2}))

figtitle = 'selectivityEmergenceDistribution';
figfile = [figDir,'PopSelectivity\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);
