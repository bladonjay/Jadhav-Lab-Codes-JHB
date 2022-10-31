%cs_numSelCellsLearning

[topDir, figDir] = cs_setPaths();
regions = {'CA1','PFC'};

figure,
hold on

pre_all = [];
post_all = [];
for r = 1:length(regions)
    region = regions{r};
    
    load([topDir,'AnalysesAcrossAnimals\selectiveCells_novel_prelearn_',region]);
    prelearn = selectivecells;
    
    load([topDir,'AnalysesAcrossAnimals\selectiveCells_novel_postlearn_',region]);
    postlearn = selectivecells;
    
    pre_all = [pre_all,size(prelearn,1)];
    post_all = [post_all,size(postlearn,1)];
    
end

bar([1:length(regions)],[pre_all;post_all]')

xticks([1:length(regions)+1]);
xticklabels(regions);
ylabel('Number of Cells');
legend({'prelearn','postlearn'},'Location','NorthWest');

figfile = [figDir,'OdorSelectivity\numSelCellsLearning'];
print('-dpdf',figfile);
print('-djpeg',figfile);