clear all

topDir = cs_setPaths();

%get only odorplace days
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


load([topDir, 'AnalysesAcrossAnimals\allCells_CA1']);
allCA1 = allcells;
load([topDir, 'AnalysesAcrossAnimals\allCells_PFC']);
allPFC = allcells;

load([topDir, 'AnalysesAcrossAnimals\runCells_CA1']);
activeCA1 = runcells;
load([topDir, 'AnalysesAcrossAnimals\runCells_PFC']);
activePFC = runcells;

load([topDir, 'AnalysesAcrossAnimals\pyrCells_CA1']);
pyrCA1 = pyrcells;
pyrCA1 = activeCA1(ismember(activeCA1,pyrCA1,'rows'),:);
load([topDir, 'AnalysesAcrossAnimals\pyrCells_PFC']);
pyrPFC = pyrcells;
pyrPFC = activePFC(ismember(activePFC,pyrPFC,'rows'),:);


load([topDir, 'AnalysesAcrossAnimals\interneurons_CA1']);
intCA1 = interneurons;
intCA1 = activeCA1(ismember(activeCA1,intCA1,'rows'),:);
load([topDir, 'AnalysesAcrossAnimals\interneurons_PFC']);
intPFC = interneurons;
intPFC = activePFC(ismember(activePFC,intPFC,'rows'),:);


load([topDir, 'AnalysesAcrossAnimals\npCells_CA1']);
npCA1 = npCells;
load([topDir, 'AnalysesAcrossAnimals\npCells_PFC']);
npPFC = npCells;

load([topDir, 'AnalysesAcrossAnimals\selectiveCells_CA1']);
selCA1 = selectivecells;
load([topDir, 'AnalysesAcrossAnimals\selectiveCells_PFC']);
selPFC = selectivecells;


sessions = unique([allCA1(:,[1,2]);allPFC(:,[1,2])],'rows');
sessions = sessions(ismember(sessions,goodsessions,'rows'),:);

%header labels: animal, day, #CA1 total, #CA1 active, #CA1 pyr, #CA1 int, #CA1 odorresp, #CA1
%selective, #PFC total, #PFC active, #PFC pyr, #PFC int, #PFC odorresp, #PFC
%selective
masterTable = zeros(size(sessions,1),14);
for s = 1:size(sessions,1)
    session = sessions(s,:);
    masterTable(s,[1,2]) = session;
    
    %CA1
    N = ismember(allCA1(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,3) = N;
    
    N = ismember(activeCA1(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,4) = N;
    
    N = ismember(pyrCA1(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,5) = N;
    
    N = ismember(intCA1(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,6) = N;
    
    N = ismember(npCA1(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,7) = N;
    
    N = ismember(selCA1(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,8) = N;
    
    %PFC
    N = ismember(allPFC(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,9) = N;
    
    N = ismember(activePFC(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,10) = N;
    
    N = ismember(pyrPFC(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,11) = N;
    
    N = ismember(intPFC(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,12) = N;
    
    N = ismember(npPFC(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,13) = N;
    
    N = ismember(selPFC(:,[1,2]),session,'rows');
    N = sum(N);
    masterTable(s,14) = N;
end

totals = sum(masterTable(:,3:end));
persession = totals./size(sessions,1);

pad = [nan(2,2),[totals;persession]];
masterTable = [masterTable;pad];

%header labels: animal, day, #CA1 total, #CA1 pyr, #CA1 int, #CA1 odorresp, #CA1
%selective, #PFC total, #PFC pyr, #PFC int, #PFC odorresp, #PFC
%selective
headings = {'Animal','Session','Total CA1','Active CA1','CA1 pyr','CA1 int','CA1 odor resp.','CA1 odor select.','Total PFC','Active PFC','PFC pyr','PFC int','PFC odor resp.','PFC odor select.'};
masterTable = [headings; num2cell(masterTable)];

xlswrite([topDir,'CellMetaData2.xlsx'], masterTable);

