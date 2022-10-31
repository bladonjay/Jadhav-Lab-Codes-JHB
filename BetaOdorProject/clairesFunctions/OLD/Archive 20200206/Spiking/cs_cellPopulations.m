%cs_cellPopulations

%Get number of cells:
%     1. all cells
%     2. cells that spike during run
%     3. cells that spike during NP
%     4. cells that are odor selective
%     5. cells with fields on track
clear
topDir = cs_setPaths();
regions = {'CA1','PFC'};

for r = 1:length(regions)
    region = regions{r};
        
    %start with all cells
        load([topDir,'AnalysesAcrossAnimals\allCells_',region]);
        cellpops.allCells = allcells;
        
        % run vs sleep
        load([topDir,'AnalysesAcrossAnimals\runCells_',region]);
        cellpops.runCells = runcells;
        
        sleepCells = allcells(~ismember(allcells,runcells,'rows'),:);
        cellpops.sleepCells = sleepCells;
        
                %test = npCells(~ismember(npCells,test,'rows'),:);

        
        %pyr vs int
        load([topDir,'AnalysesAcrossAnimals\pyrCells_',region]);
        pyrcells = runcells(ismember(runcells,pyrcells,'rows'),:); %run cells that are pyr
        cellpops.pyrCells = pyrcells;
        intCells = runcells(~ismember(runcells,pyrcells,'rows'),:);
        cellpops.intCells = intCells;
        
        
        
        % np vs non-np
        load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
        cellpops.npCells = npCells;
        npCells = pyrcells(ismember(pyrcells,npCells,'rows'),:);
        nonNPCells = pyrcells(~ismember(pyrcells,npCells,'rows'),:);
        cellpops.nonNPCells = nonNPCells;
        
        
        %selective vs nonselective
        load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region]);
        cellpops.selectiveCells = selectivecells;
        nonselective = npCells(~ismember(npCells,selectivecells,'rows'),:);
        cellpops.nonselectiveCells = nonselective;

        save([topDir,'AnalysesAcrossAnimals\cellPopulations_',region],'cellpops');
        clear cellpops
end