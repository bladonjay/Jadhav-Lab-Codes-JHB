%identifies cells that spiked during rewards. Excludes ripples. 
clear
animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC'};
topDir = cs_setPaths;

for r = 1:length(regions)
    region = regions{r};
    
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal,'Expt\',animal,'_direct\'];
        load([animDir, animal, 'tetinfo.mat']);
        files = dir([animDir,animal,'rewards*']);
        
        for d = 1:length(files)
            load([animDir,files(d).name]);
            day = find(~cellfun(@isempty, rewards));
            epochs = find(~cellfun(@isempty, rewards{day}));
            daystr = getTwoDigitNumber(day);
            
            load([animDir, animal, 'ripples',daystr,'.mat']);
            load([animDir, animal, 'spikes',daystr,'.mat']);
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                leftRewards = rewards{day}{epoch}.leftWindows;
                rightRewards = rewards{day}{epoch}.rightWindows;
                
                %find tetrode with most cells, use as riptet
                tetfilter = '(isequal($area,''CA1''))';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                
                numcells = [];
                for t = 1:length(tets)
                    if ~isempty(tetinfo{day}{epoch}{t})
                        num = tetinfo{day}{epoch}{t}.numcells;
                    else
                        num = 0;
                    end
                    numcells = [numcells;num];    
                end
                [~,ind] = max(numcells);
                tet = tets(ind);
                
                %epcells = 
                %ripTimes = 
            end
        end
    end
    %save in AnalysesAcrossAnimals
end