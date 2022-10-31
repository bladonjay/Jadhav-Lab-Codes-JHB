clear 
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
%animals = {'CS44'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];
learningtypes = {'prelearn','postlearn'};


for r = 1:length(regions)
    region = regions{r};
    
    
    
    for l = 1:length(learningtypes)
        learning = learningtypes{l};
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
        
        days1 = cs_getRunEpochs(animDir,animal,'novelodor');
        days2 = cs_getRunEpochs(animDir,animal,'novelodor2');
        dayeps = unique([days1; days2],'rows');
        days = unique(dayeps(:,1));
        
        
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'nosepokeWindow',daystr,'.mat'])
            load([animDir, animal,'odorTriggers',daystr,'.mat']);
            
            runeps = dayeps((dayeps(:,1) == day),2);
            %runeps = find(~cellfun(@isempty,nosepokeWindow{day}));
            
           
            for c = 1:size(daycells,1)

                npspikes = 0;
                totaltrigs = 0;
                cell = daycells(c,:);
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    
                    if ~isfield(odorTriggers{day}{epoch},learning)
                        continue
                    end
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data) 
                            epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                            
                            [cl, cr, il, ir] = cs_getSpecificTrialTypeInds_learning(odorTriggers{day}{epoch},learning);
                            
                            trigs = nosepokeWindow{day}{epoch}([cl;cr;il;ir],:);
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
    
    save([dataDir,'npCells_novel_',learning,'_',region,'.mat'], 'npCells')
    end
    end