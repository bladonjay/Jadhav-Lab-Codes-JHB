%Overlap
[topDir, figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
figDir = [figDir,'NicePPTFigures\'];
freqbands = {'beta','theta'};
eegregions = {'CA1','PFC','OB'};

for f = 1:length(freqbands)
    freqband = freqbands{f};
    
    for r = 1:length(eegregions)
        eegregion = eegregions{r};
 
        region = 'CA1';
        
        load([dataDir, 'npCells_',region,'.mat']);
        load([dataDir, 'selectiveCells_',region,'.mat']);
        load([dataDir, 'plCells_',freqband,'_',region,'-',eegregion,'.mat']);
        
        
        bothCA1 = intersect(selectivecells,plcells,'rows');
        
        fractionselective = size(selectivecells,1) / size(npCells,1);
        fractionpl = size(plcells,1)/size(npCells,1);
        fractionboth = size(bothCA1,1)/size(npCells,1);
        
        probsel = binofit(size(selectivecells,1),size(npCells,1),0.5);
        probpl = binofit(size(plcells,1),size(npCells,1),0.5);
        probbothCA1 = binofit(size(bothCA1,1),size(npCells,1),(probpl*probsel));
        
        
        barsCA1 = [fractionselective,fractionpl,fractionboth];
        barsCA1 = barsCA1*100; %percentage
        
        region = 'PFC';
        
        load([dataDir, 'npCells_',region,'.mat']);
        load([dataDir, 'selectiveCells_',region,'.mat']);
        load([dataDir, 'plCells_',freqband,'_',region,'-',eegregion,'.mat']);
        
        bothPFC = intersect(selectivecells,plcells,'rows');
        
        fractionselective = size(selectivecells,1) / size(npCells,1);
        fractionpl = size(plcells,1)/size(npCells,1);
        fractionboth = size(bothPFC,1)/size(npCells,1);
        
        barsPFC = [fractionselective,fractionpl,fractionboth];
        barsPFC = barsPFC*100; %percentage
        
        probsel = binofit(size(selectivecells,1),size(npCells,1),0.5);
        probpl = binofit(size(plcells,1),size(npCells,1),0.5);
        probbothPFC = binofit(size(bothPFC,1),size(npCells,1),(probpl*probsel));
        
        
        figure,
        bar([1,2],[barsCA1;barsPFC])
        box off
        ylabel('Percentage of Cells');
        xticklabels({'CA1','PFC'})
        axis([0.5 2.5, 0 100])
       
        figfile = [figDir,'SelectivityPLOverlap_',eegregion,'_', freqband];
        
        legend({'Selective Cells',[eegregion, ' ',freqband, ' Phase Locked Cells'],'Overlap'});
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        saveas(gcf,figfile,'fig');
    end
end
