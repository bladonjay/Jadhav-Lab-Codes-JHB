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
        
     load([topDir,'AnalysesAcrossAnimals\runCells_',region]);
        allcells = runcells;
        cellpops.runCells = allcells;
        
   
         sessions = unique(allcells(:,[1,2]),'rows');
         sessions = sessions(ismember(sessions,goodsessions,'rows'),:);     
        
        % int
        load([topDir,'AnalysesAcrossAnimals\pyrCells_',region]);
%         N = ismember(interneurons(:,[1,2]),sessions,'rows');
%         interneurons = interneurons(N,:);
        interneurons = runcells(~ismember(runcells,pyrcells,'rows'),:);
        cellpops.intCells = interneurons;
%         
%     %start with all ints
%         load([topDir,'AnalysesAcrossAnimals\interneurons_',region]);
%         
%         sessions = unique(interneurons(:,[1,2]),'rows');
%         sessions = sessions(ismember(sessions,goodsessions,'rows'),:);
%         N = ismember(interneurons(:,[1,2]),sessions,'rows');
%         interneurons = interneurons(N,:);
%         cellpops.allints = interneurons;
        
%         %pyr vs int
%         load([topDir,'AnalysesAcrossAnimals\pyrCells_',region]);
%         N = ismember(pyrcells(:,[1,2]),sessions,'rows');
%         pyrcells = pyrcells(N,:);
%         cellpops.pyrCells = pyrcells;
%         
%         intCells = allcells(~ismember(allcells,pyrcells,'rows'),:);
%         N = ismember(intCells(:,[1,2]),sessions,'rows');
%         intCells = intCells(N,:);
%         cellpops.intCells = intCells;
        
%         % run vs sleep
%         load([topDir,'AnalysesAcrossAnimals\runCells_',region]);
%         runcells = pyrcells(ismember(pyrcells,runcells,'rows'),:);
%         cellpops.runCells = runcells;
%         
%         sleepCells = pyrcells(~ismember(pyrcells,runcells,'rows'),:);
%         cellpops.sleepCells = sleepCells;
        
        
        % np vs non-np
        load([topDir,'AnalysesAcrossAnimals\npInt_',region]);
        cellpops.npCells = npInt;
        npInt = interneurons(ismember(interneurons,npInt,'rows'),:);
        nonNPCells = interneurons(~ismember(interneurons,npInt,'rows'),:);
        cellpops.nonNPCells = nonNPCells;
        
        
        %selective vs nonselective
        load([topDir,'AnalysesAcrossAnimals\selectiveINs_',region]);
        cellpops.selectiveCells = selectiveint;
        nonselective = npInt(~ismember(npInt,selectiveint,'rows'),:);
        cellpops.nonselectiveCells = nonselective;

        save([topDir,'AnalysesAcrossAnimals\intPopulations_',region],'cellpops');
        clear cellpops
end