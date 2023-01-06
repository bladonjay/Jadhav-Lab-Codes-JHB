
% Saving and Plotting Variables

% % 1) PCA_AvgDiscrT
% % 8X2, 8 animals, one column each for CA1 and PFC
% 
% load('PCA_AvgDiscrT.mat');
% figure; hold on; boxplot(PCA_AvgDiscrT)
% %x=repmat(1:2,8,1); 
% scatter(x(:),PCA_AvgDiscrT(:),'filled','MarkerFaceAlpha',1','jitter','on','jitterAmount',0.15);
% 
% nanmean(PCA_AvgDiscrT)
% nanmedian(PCA_AvgDiscrT)
%p = kruskalwallis(PCA_AvgDiscrT)




% 2) PCA_AvgDiscrT_allsess
% 39 rows, one for each day
% Columns:
% 1: Anim
% 2: Sess
% 3: CA1
% 4: PFC
% Will have NaNs (Total 24 non-NaN sessions)
try
    load('E:\Brandeis datasets\OdorPlaceAssociation\AnalysesAcrossAnimals\PCA_AvgDiscrT_allsess.mat');
catch
    [myfile,mydir]=uigetfile('find the you created when you ran ''Run_PCA_SJ_copy''');
    load(fullfile(mydir,myfile));
end
CA1_discr = PCA_AvgDiscrT_allsess(~isnan(PCA_AvgDiscrT_allsess(:,3)),3); % 24 sessions/ ori26
PFC_discr = PCA_AvgDiscrT_allsess(~isnan(PCA_AvgDiscrT_allsess(:,4)),4); % 24 sessions/ ori25

x1=ones(length(CA1_discr),1); 
x2=2*ones(length(PFC_discr),1); 
figure; hold on; 
%boxplot([CA1_discr,PFC_discr]);
bar(1,mean(PCA_AvgDiscrT_allsess(:,3),'omitnan'),.6,'LineStyle','none','FaceAlpha',.5,'FaceColor',rgbcolormap('tomato'));
hold on;
bar(2,mean(PCA_AvgDiscrT_allsess(:,4),'omitnan'),.6,'LineStyle','none','FaceAlpha',.5,'FaceColor',rgbcolormap('aqua'));

hold on;
%plot(PCA_AvgDiscrT_allsess(:,3:4)','b','LineWidth',2)
errorbar([1 2],mean(PCA_AvgDiscrT_allsess(:,[3 4]),'omitnan'),nanstd(PCA_AvgDiscrT_allsess(:,[3 4])),'.','CapSize',0,'LineWidth',3,...
    'Color','k')

scatter(x1,CA1_discr,24,rgbcolormap('tomato'),'filled','MarkerFaceAlpha',1,'jitter','on','jitterAmount',0.15);
scatter(x2,PFC_discr,24,rgbcolormap('teal'),'filled','MarkerFaceAlpha',1,'jitter','on','jitterAmount',0.15);
ylabel('PC Discrimination time (mSec)'); set(gca,'XTick',[1 2],'XTickLabel',{'CA1','PFC'});
% Need to make a paired plot
%--------------------------

mean(CA1_discr), mean(PFC_discr)
median(CA1_discr), median(PFC_discr)
title(sprintf('Signed-rank P= %.2e',signrank(CA1_discr,PFC_discr)));


% 3) Example CA1 Discrimination. 3 PCs, 12 cells

load PCA_CA1Popln_example1.mat 
figure; hold on;
title ('PCA Distance Metric, 3 PCs; CA1 example');
plot(timeaxis,store_avg_DISTPC,'ko-', 'LineWidth',4,'MarkerSize',16);
plot(timeaxis,shuffavgDIST_PC,'ko-', 'LineWidth',8);
plot(timeaxis,shuffavgCIDIST_PC(2,:),'k--', 'LineWidth',2); % 95% CI
plot(timeaxis,shuffavgCIDIST_PC(1,:),'k--', 'LineWidth',2); % 95% CI 
figure;
plot(timeaxis,store_avg_DISTPC,'ko-', 'LineWidth',4,'MarkerSize',16);


% 4) Example PFC Discrimination. 3 PCs, 9 cells

load PCA_PFCPopln_example1.mat 
figure; hold on;
title ('PCA Distance Metric, 3 PCs; PFC example');
plot(timeaxis,store_avg_DISTPC,'ko-', 'LineWidth',4,'MarkerSize',16);
plot(timeaxis,shuffavgDIST_PC,'ko-', 'LineWidth',8);
plot(timeaxis,shuffavgCIDIST_PC(2,:),'k--', 'LineWidth',2); % 95% CI
plot(timeaxis,shuffavgCIDIST_PC(1,:),'k--', 'LineWidth',2); % 95% CI 




% 5) Example: Raw PCA left vs PCA right 3d example - from which distance is
% computed

load PCA_raweg2.mat
figure; hold on;
plot3(meanleft(:,1),meanleft(:,2),meanleft(:,3),'ro-', 'LineWidth',4);
plot3(meanright(:,1),meanright(:,2),meanright(:,3),'bx-', 'LineWidth',4);




% 6) Example: Single trial PCA vector distance and corresponding latency

load single_trial_PCA_latency_eg1.mat
figure; hold on; 

% Plot Avg
plot(timeaxis,store_avg_DISTPC,'ko-', 'LineWidth',4,'MarkerSize',16);
plot(timeaxis,shuffavgDIST_PC,'ko-', 'LineWidth',8);
plot(timeaxis,shuffavgCIDIST_PC(2,:),'k--', 'LineWidth',2); % 95% CI
plot(timeaxis,shuffavgCIDIST_PC(1,:),'k--', 'LineWidth',2); % 5% CI

% Plot Early trial
plot(timeaxis,currdis,'ro-', 'LineWidth',2,'MarkerSize',8);
plot(currlat,0,'rsq','MarkerSize',16);
title ('Single Trial PC Distance Metric and Dec. Latency ');

% Plot Late trial
plot(timeaxis,currdis2,'bo-', 'LineWidth',2,'MarkerSize',8);
plot(currlat2,0,'bsq','MarkerSize',16);

