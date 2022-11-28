%make list of selective cells

%topDir = 'F:\Data\OdorPlaceAssociation\';
%topDir = 'D:\OdorPlaceAssociation\';
clear
topDir = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
dataDir = [topDir,'AnalysesAcrossAnimals\'];

for r = 1:length(regions)
region = regions{r};
    selectiveint = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['((isequal($area,''',region, ... 
            ''')) && strcmp($type,''int'') && ',...
            '(strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
        
        cells = evaluatefilter(cellinfo,cellfilter);
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        animvector = repmat(a, size(cells,1), 1);
        cells = [animvector, cells];
        selectiveint = [selectiveint; cells];
    end
    fprintf('there are %d selective INs in %s \n', size(selectiveint,1), regions{r});
    save([dataDir,'selectiveINs_',region,'.mat'], 'selectiveint');
end


        