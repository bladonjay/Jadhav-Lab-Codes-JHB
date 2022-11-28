%make list of place cells

clear
topDir = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
dataDir = [topDir,'AnalysesAcrossAnimals\'];

for r = 1:length(regions)
region = regions{r};
    placecells = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['isequal($area,''',region,''') && strcmp($placetag, ''placecell'')'];
        cells = evaluatefilter(cellinfo,cellfilter);
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        animvector = repmat(a, size(cells,1), 1);
        cells = [animvector, cells];
        placecells = [placecells; cells];
    end
    
    save([dataDir,'placeCells_',region,'.mat'], 'placecells');
end


        