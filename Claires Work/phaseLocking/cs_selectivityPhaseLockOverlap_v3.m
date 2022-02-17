%Overlap
[topDir, figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
figDir = [figDir,'PhaseLocking\'];
freqbands = {'beta'};
eegregions = {'CA1','PFC'};

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
        load([dataDir, 'PhaseLocking\plCells_',freqband,'_',region,'-',eegregion,'.mat']);
        allplcells = [allplcells; plcells];
    end
    plcells = unique(allplcells,'rows');
    
    bothCA1 = intersect(selectivecells,plcells,'rows');
    
    fractionselective = size(selectivecells,1) / size(npCells,1);
    fractionpl = size(plcells,1)/size(npCells,1);
    fractionboth = size(bothCA1,1)/size(npCells,1);
    
    probbothCA1 = (fractionpl*fractionselective)*100;
    %probbothCA1 = fractionselective*fractionpl*100;
    pboth_CA1 = binopdf(size(bothCA1,1),size(npCells,1),(fractionpl*fractionselective));
    
    barsCA1 = [fractionselective,fractionpl,fractionboth];
    barsCA1 = barsCA1*100; %percentage
    
    X = binoinv([0.05 0.95],size(npCells,1),(fractionpl*fractionselective));
    CI_CA1 = X/size(npCells,1)*100;
   
    
    
    %% ---
    region = 'PFC';
    
    load([dataDir, 'npCells_',region,'.mat']);
    load([dataDir, 'selectiveCells_',region,'.mat']);
    %         load([dataDir, 'plCells_',freqband,'_',region,'-',eegregion,'.mat']);
    
    %count cells as phase locked as long as they are phase locked to
    %beta from at least one region
    allplcells = [];
    for r = 1:length(eegregions)
        eegregion = eegregions{r};
        load([dataDir, 'PhaseLocking\plCells_',freqband,'_',region,'-',eegregion,'.mat']);
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
    probbothPFC = (fractionpl*fractionselective)*100;
    %probbothCA1 = fractionselective*fractionpl*100;
    pboth_PFC = binopdf(size(bothPFC,1),size(npCells,1),(fractionpl*fractionselective));
    %probbothPFC = fractionselective*fractionpl*100;
    X = binoinv([0.05 0.95],size(npCells,1),(fractionpl*fractionselective));
    CI_PFC = X/size(npCells,1)*100;
    
    disp(['CA1 p =',num2str(pboth_CA1)]);
    disp(['PFC p =',num2str(pboth_PFC)]);

    
    figure,
    bar([1,2],[barsCA1;barsPFC])
    box off
    ylabel('Percentage of Cells');
    xticklabels({'CA1','PFC'})
    axis([0.5 2.5, 0 100])
    hold on
    plot([1.1 1.35],[probbothCA1 probbothCA1], 'k--')
    plot([2.1 2.35],[probbothPFC probbothPFC], 'k--')
    plot([1.25 1.25],CI_CA1, 'r--')
    plot([2.25 2.25],CI_PFC, 'r--');
    text(1,70,{num2str(barsCA1(1));num2str(barsCA1(2));num2str(barsCA1(3));num2str(barsPFC(1));num2str(barsPFC(2));num2str(barsPFC(3))})
    
    figfile = [figDir,'SelectivityPLOverlap_',freqband];
    
    legend({'Selective Cells',[freqband, ' Phase Locked Cells'],'Overlap'});
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
    saveas(gcf,figfile,'fig');
    %     end
end
