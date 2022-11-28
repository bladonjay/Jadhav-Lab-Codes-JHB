%Only Looking at cells that phase lock just before the turn away from NP
clear

[topDir, figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
%figDir = [figDir,'NicePPTFigures\'];
freqbands = {'beta'};
eegregions = {'CA1','PFC','OB'};

for f = 1:length(freqbands)
    freqband = freqbands{f};
    
    region = 'CA1';
    
    load([dataDir, 'npCells_',region,'.mat']);
    load([dataDir, 'selectiveCells_',region,'.mat']);
    %load([dataDir, 'plCells_',freqband,'_',region,'-',eegregion,'.mat']);
    
    %count cells as phase locked as long as they are phase locked to
    %beta from at least one region
    allplcells = [];
    for r = 1:length(eegregions)
        eegregion = eegregions{r};
        load([dataDir,'PhaseLocking\NPChoice\',region,'-',eegregion,'.mat']);
        plcells = unique(plcells(:,[1 2 4 5]),'rows');
        allplcells = [allplcells; plcells];
    end
    plcells = unique(allplcells,'rows');
    
    bothCA1 = intersect(selectivecells,plcells,'rows');
    
    fractionselective = size(selectivecells,1) / size(npCells,1);
    fractionpl = size(plcells,1)/size(npCells,1);
    fractionboth = size(bothCA1,1)/size(npCells,1);
    
%     probsel = binofit(size(selectivecells,1),size(npCells,1),0.5);
%     probpl = binofit(size(plcells,1),size(npCells,1),0.5);
%     probbothCA1 = binofit(size(bothCA1,1),size(npCells,1),(probpl*probsel));
    probbothCA1 = fractionselective*fractionpl*100;
    
    barsCA1 = [fractionselective,fractionpl,fractionboth];
    barsCA1 = barsCA1*100; %percentage
    
    region = 'PFC';
    
    load([dataDir, 'npCells_',region,'.mat']);
    load([dataDir, 'selectiveCells_',region,'.mat']);
    %         load([dataDir, 'plCells_',freqband,'_',region,'-',eegregion,'.mat']);
    
    %count cells as phase locked as long as they are phase locked to
    %beta from at least one region
    allplcells = [];
    for r = 1:length(eegregions)
        eegregion = eegregions{r};
        load([dataDir,'PhaseLocking\NPChoice\',region,'-',eegregion,'.mat']);
        plcells = unique(plcells(:,[1 2 4 5]),'rows');
        allplcells = [allplcells; plcells];
    end
    plcells = unique(allplcells,'rows');
    
    bothPFC = intersect(selectivecells,plcells,'rows');
    
    fractionselective = size(selectivecells,1) / size(npCells,1);
    fractionpl = size(plcells,1)/size(npCells,1);
    fractionboth = size(bothPFC,1)/size(npCells,1);
    
    barsPFC = [fractionselective,fractionpl,fractionboth];
    barsPFC = barsPFC*100; %percentage
    
%     probsel = binofit(size(selectivecells,1),size(npCells,1),0.5);
%     probpl = binofit(size(plcells,1),size(npCells,1),0.5);
%     probbothPFC = binofit(size(bothPFC,1),size(npCells,1),(probpl*probsel));

    probbothPFC = fractionselective*fractionpl*100;
    
    
    
    figure,
    bar([1,2],[barsCA1;barsPFC])
    box off
    ylabel('Percentage of Cells');
    xticklabels({'CA1','PFC'})
    axis([0.5 2.5, 0 100])
    hold on
    plot([1.1 1.35],[probbothCA1 probbothCA1], 'k--')
    plot([2.1 2.35],[probbothPFC probbothPFC], 'k--')
    
    text(1,70,{num2str(barsCA1(1));num2str(barsCA1(2));num2str(barsCA1(3));num2str(barsPFC(1));num2str(barsPFC(2));num2str(barsPFC(3))})
    
    %figfile = [figDir,'SelectivityPLOverlap_',freqband,'_NPChoice'];
    
    legend({'Selective Cells',[freqband, ' Phase Locked Cells'],'Overlap'});
    
    %print('-djpeg', figfile);
    %print('-dpdf', figfile);
    %saveas(gcf,figfile,'fig');
    %     end
end
