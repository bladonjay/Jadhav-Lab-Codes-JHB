clear 
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
%regions = {'CA1'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
    region = regions{r};
    
    npCells =[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['isequal($area,''',region,''') && isequal($type, ''pyr'')'];
%         cellfilter = ['isequal($type, ''pyr'')']
%         cellfilter = ['isequal($area,''',region,''')']
        cells = evaluatefilter(cellinfo,cellfilter);
        
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        runEps = cs_getRunEpochs(animDir,animal,'odorplace');
        days = unique(runEps(:,1));
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'nosepokeWindow',daystr,'.mat'])
            eps = runEps((runEps(:,1) == day),2);
            
           
            for c = 1:size(daycells,1)

                npspikes = 0;
                totaltrigs = 0;
                cell = daycells(c,:);
                
                for ep = 1:length(eps)
                    epoch = eps(ep);
                    
                    
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data) 
                            epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);

                            trigs = nosepokeWindow{day}{epoch};
                            totaltrigs = totaltrigs + size(trigs,1);
                            
                            winspikes = epspikes(isExcluded(epspikes, trigs));
                            npspikes = npspikes + length(winspikes);
                           
                        end

                        
                        
%                         for t = 1:size(trigs,1)
%                         trigwin = trigs(t,:);
%                         winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
%                         npspikes = npspikes + length(winspikes);
% 
%                         end
                        
                    end

                end
                
                avgSpikesPerTrial = npspikes/totaltrigs;
                
                if avgSpikesPerTrial >= 1
                   npCells = [npCells; a, cell];
                end
            end
        end
    end
    
    save([dataDir,'npCells_',region,'.mat'], 'npCells')
end