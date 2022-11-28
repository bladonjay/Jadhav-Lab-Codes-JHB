function cs_plotPlSelectivityOverlap(dataDir, figDir)

%cs_plotPlSelectivityOverlap('E:\AnalysesAcrossAnimals\','E:\Figures\', {'CA1','PFC'})

%dataDir = 'E:\AnalysesAcrossAnimals\';
%figDir = 'E:\Figures\';

regions = {'CA1','PFC'};
winstr = '_-500-1000ms';

load(['phaseLockedCellFractions.mat'])

datatoplot = zeros(2,3);
for r = 1:length(regions)
    load([dataDir,'selectivityFractions_',regions{r},winstr,'.mat'])
    load([dataDir,'phaseLockedSelectiveCells_',regions{r},'.mat'])
    load([dataDir,'cellSelectivityData_',regions{r},winstr,'.mat']);
    fractselective = selectivityFractions.fractSelective;
    
    eval(['fractphaselocked = ',regions{r},'phaseLockedFraction;'])
    
    fractoverlap = length(phaseLockedSelectiveCells)/length(cellSelectivity);
    
    
    datatoplot(r,1:3) = [fractphaselocked, fractselective, fractoverlap];


end

%save([dataDir,'PL_SelectiveCells.mat'],'PL_SelectiveCells');

figure
b = bar([1,2], datatoplot);

c = [ rgb('LightSteelBlue'); rgb('SlateGrey'); rgb('Khaki')];
b(1).FaceColor = c(1,:);
b(2).FaceColor = c(2,:);
b(3).FaceColor = c(3,:);



%set(gca,'fontsize',20);


xticklabels({'CA1','PFC'})
ylabel('Fraction of Cells')

set(gca,'fontsize',32);
box off
set(gcf, 'Position', [50 50 800 800]);

l = legend({'Phase Locked Cells', 'Selective Cells', 'Selective-Phase Locked Cells'},'Position',[.16 .75 1 .1]);
l.FontSize = 25;
legend('boxoff')

figtitle = 'OverlappingCells';
figfile = [figDir,'NicePPTFigures\',figtitle];
%saveas(gcf,figfile,'fig');
print( '-dpdf', figfile);
print('-djpeg', figfile);