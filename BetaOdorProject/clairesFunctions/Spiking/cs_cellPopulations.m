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
eegregions={'CA1','PFC','OB'};

for r = 1:length(regions)
    region = regions{r};
    
    %start with all cells, at least 100 spikes
    load([topDir,'AnalysesAcrossAnimals\allCells_',region]);
    cellpops.allCells = allcells;
    
    % run vs sleep: at least 1 spike in run
    load([topDir,'AnalysesAcrossAnimals\runCells_',region]);
    cellpops.runCells = runcells;
    
    sleepCells = allcells(~ismember(allcells,runcells,'rows'),:);
    cellpops.sleepCells = sleepCells;
    
    %test = npCells(~ismember(npCells,test,'rows'),:);
    
    
    %run pyrs and run INs 
    load([topDir,'AnalysesAcrossAnimals\pyrCells_',region]);
    pyrCells = runcells(ismember(runcells,pyrcells,'rows'),:); %run cells that are pyr
    cellpops.pyrCells = pyrCells;
    intCells = runcells(~ismember(runcells,pyrCells,'rows'),:);
    cellpops.intCells = intCells;
    
    
    
    % np pyrs and non np pyrs
    %load([topDir,'AnalysesAcrossAnimals\npCells_',region,'_old']);
    load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
    npPyrams= pyrCells(ismember(pyrCells,npCells,'rows'),:);
    cellpops.npPyrams = npPyrams;
    nonNPCells = pyrCells(~ismember(pyrCells,npPyrams,'rows'),:);
    cellpops.nonNPCells = nonNPCells;
    clear npCells;
    
    % np vs non-np INs
    %load([topDir,'AnalysesAcrossAnimals\npInt_',region,'_old']);
    load([topDir,'AnalysesAcrossAnimals\npInt_',region]);
    npInt= intCells(ismember(intCells,npInt,'rows'),:);
    cellpops.npInt = npInt;
    nonNPInt = pyrCells(~ismember(intCells,npInt,'rows'),:);
    cellpops.nonNPInt = nonNPInt;
    
    
    %selective vs nonselective Pyrams 
    load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region, '_new']);
    selectivePyram= npPyrams(ismember(npPyrams,selectivecells,'rows'),:);
    cellpops.selectiveCells = selectivePyram; % NOT TRUE, these are pyrams and INTs
    nonselective = npPyrams(~ismember(npPyrams,selectivecells,'rows'),:);
    cellpops.nonselectiveCells = nonselective;
    
     %selective vs nonselective Int 
    load([topDir,'AnalysesAcrossAnimals\selectiveINs_',region]);
    selectiveINs= npInt(ismember(npInt,selectivecells,'rows'),:);
    cellpops.selectiveCells = selectivePyram; % NOT TRUE, these are pyrams and INTs
    nonselective = npPyrams(~ismember(npPyrams,selectivecells,'rows'),:);
    cellpops.nonselectiveCells = nonselective;
    
    
    
    
    % now beta phaselocked pyrams from npCells
    % use npPyrams and npInt
    % have to pull data from all three and then can union them
    
    plAll=[];
    for er=1:length(eegregions)
        load([topDir,'AnalysesAcrossAnimals\PhaseLocking\plCells_',...
            'beta_',regions{r},'-',eegregions{er},'.mat'], 'plcells');
        plAll=[plAll; plcells];
        cellpops.(['betaPhaseLockedTo', eegregions{er}])=plcells;
    end
       plAll=unique(plAll,'rows');
       cellpops.betaPhaseLockedAny=plAll;

    
    plAll=[];
    for er=1:length(eegregions)
        load([topDir,'AnalysesAcrossAnimals\PhaseLocking\plCells_',...
            'resp_',regions{r},'-',eegregions{er},'.mat'], 'plcells');
        plAll=[plAll; plcells];
        cellpops.(['respPhaseLockedTo', eegregions{er}])=plcells;
    end
       plAll=unique(plAll,'rows');
       cellpops.respPhaseLockedAny=plAll;
    
    
    save([topDir,'AnalysesAcrossAnimals\cellPopulations_',region],'cellpops');
    clear cellpops
end
%%
clear
topDir = cs_setPaths();
regions = {'CA1','PFC'};

for r=1:length(regions)
    load([topDir,'AnalysesAcrossAnimals\cellPopulations_',regions{r}],'cellpops');
    eval(['cellpops_',regions{r}, '=cellpops; clear cellpops']);
end
openvar('cellpops_CA1')