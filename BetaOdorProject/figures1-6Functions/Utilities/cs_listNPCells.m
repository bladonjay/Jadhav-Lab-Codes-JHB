
win = [0 1];
animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
region = regions{r};
    
    npCells =[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['isequal($area,''',region,''')'];
        
        cells = evaluatefilter(cellinfo,cellfilter);
        
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        days = unique(cells(:,1));
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
           
            for c = 1:size(daycells,1)

                npspikes = 0;
                
                cell = daycells(c,:);
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data) 
                            epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                        end

                        trigs = odorTriggers{day}{epoch}.allTriggers;
                        
                        for t = 1:length(trigs)
                        trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
                        winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                        npspikes = npspikes + length(winspikes);

                        end
                        
                    end

                end
                
                if npspikes >= length(trigs)
                   npCells = [npCells; a, cell];
                end
            end
        end
    end
    
    save([dataDir,'npCells_',region,'.mat'], 'npCells')
end