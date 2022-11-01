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

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
goodsessions = [];
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    days = unique(runeps(:,1));
    an = repmat(a,length(days),1);
    goodsessions = [goodsessions;an,days];
end

for r = 1:length(regions)
    region = regions{r};

    % runCells are those that spike at least once during run epochs
    % uses function CS_listrRunCells
    load([topDir,'AnalysesAcrossAnimals\runCells_',region]);
    allcells = runcells;
    % has to spike at least once during run epochs
    cellpops.runCells = allcells;

    sessions = unique(allcells(:,[1,2]),'rows');
    sessions = sessions(ismember(sessions,goodsessions,'rows'),:);


    %pyr vs int, using the filter function (
    load([topDir,'AnalysesAcrossAnimals\pyrCells_',region]);

    pyrcells = allcells(ismember(allcells,pyrcells,'rows'),:);
    cellpops.pyrCells = pyrcells;

    intCells = allcells(~ismember(allcells,pyrcells,'rows'),:);
    N = ismember(intCells(:,[1,2]),sessions,'rows');
    intCells = intCells(N,:);
    intCells = runcells(ismember(runcells,intCells,'rows'),:);
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