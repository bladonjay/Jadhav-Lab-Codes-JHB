%make list of selective cells

%topDir = 'F:\Data\OdorPlaceAssociation\';
%topDir = 'D:\OdorPlaceAssociation\';
clear
topDir = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
dataDir = [topDir,'AnalysesAcrossAnimals\'];

for r = 2:length(regions)
region = regions{r};
    selectivecells = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['((isequal($area,''',region, ... 
            ''')) && strcmp($tag,''accepted'') && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
        
        cells = evaluatefilter(cellinfo,cellfilter);
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        animvector = repmat(a, size(cells,1), 1);
        cells = [animvector, cells];
        selectivecells = [selectivecells; cells];
    end
    
    save([dataDir,'selectiveCells_',region,'.mat'], 'selectivecells');
end


        