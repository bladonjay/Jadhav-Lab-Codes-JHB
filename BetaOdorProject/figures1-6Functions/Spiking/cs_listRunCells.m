%cs_listRunSpikes
% has to spike at least once during the odor period
win = [0 1];
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
    region = regions{r};
    
    runcells =[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        runEps = cs_getRunEpochs(animDir, animal,'odorplace');
        
        days = unique(runEps(:,1));
        
        for day = days'
            eps = runEps((runEps(:,1) == day),2);
            daycells = [];
            for ep = eps'
                
                cellfilter = ['isequal($area,''',region,''') & ($numspikes > 1)'];
                epcells = evaluatefilter(cellinfo{day}{ep},cellfilter);
                
                %noeps = animcells(:,[1 3 4]);
                %animcells = unique(noeps,'rows');
                
                %days = unique(animcells(:,1));
                
                animvect = repmat([a, day], size(epcells,1), 1);
                epcells = [animvect, epcells];
                daycells = [daycells;epcells];
                
            end
            daycells = unique(daycells,'rows');
            runcells = [runcells; daycells];
        end
    end
    
    save([dataDir,'runCells_',region,'.mat'], 'runcells')
end