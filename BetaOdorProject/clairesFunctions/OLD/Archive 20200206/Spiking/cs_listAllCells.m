
win = [0 1];
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
region = regions{r};
    
    allcells =[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['isequal($area,''',region,''') & ($numspikes > 100)']; 
        %cellfilter = ['isequal($area,''',region,''') & ($numspikes > 100)']; 
        animcells = evaluatefilter(cellinfo,cellfilter);
        
        noeps = animcells(:,[1 3 4]);
        animcells = unique(noeps,'rows');
        
        days = unique(animcells(:,1));
       
        animvect = repmat(a, size(animcells,1), 1);
        animcells = [animvect, animcells];
        
        allcells = [allcells; animcells];
    end
    
    save([dataDir,'allCells_',region,'.mat'], 'allcells')
end