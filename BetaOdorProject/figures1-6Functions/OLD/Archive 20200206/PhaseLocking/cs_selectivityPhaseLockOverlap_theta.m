%Overlap

[topDir, figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
region = 'CA1';
eegregion = 'OB';

load([dataDir, 'npCells_',region,'.mat']);
load([dataDir, 'selectiveCells_',region,'.mat']);
load([dataDir, 'plCells_theta_',region,'-',eegregion,'.mat']);

bothCA1 = intersect(selectivecells,plcells,'rows');

fractionselective = size(selectivecells,1) / size(npCells,1);
fractionpl = size(plcells,1)/size(npCells,1);
fractionboth = size(bothCA1,1)/size(npCells,1);

barsCA1 = [fractionselective,fractionpl,fractionboth];
barsCA1 = barsCA1*100; %percentage

region = 'PFC';

load([dataDir, 'npCells_',region,'.mat']);
load([dataDir, 'selectiveCells_',region,'.mat']);
load([dataDir, 'plCells_theta_',region,'-',eegregion,'.mat']);

bothPFC = intersect(selectivecells,plcells,'rows');

fractionselective = size(selectivecells,1) / size(npCells,1);
fractionpl = size(plcells,1)/size(npCells,1);
fractionboth = size(bothPFC,1)/size(npCells,1);

barsPFC = [fractionselective,fractionpl,fractionboth];
barsPFC = barsPFC*100; %percentage


figure,
bar([1,2],[barsCA1;barsPFC])
box off
ylabel('Percentage of Cells');
xticklabels({'CA1','PFC'})

figDir = [figDir,'NicePPTFigures\'];
figfile = [figDir,'SelectivityPLOverlap_RR'];

legend({'Selective Cells','Phase Locked Cells','Overlap'});
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
    saveas(gcf,figfile,'fig');
