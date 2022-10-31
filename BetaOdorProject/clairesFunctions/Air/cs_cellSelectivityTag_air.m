
%combine spiketimes across run epochs
%combine trig times across run epochs

%tag cells with selectivity index
%compare to shuffled distribution to get significance

clear
topDir = cs_setPaths;

animals = {'CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
iterations = 1000;


for r = 1:length(regions)
    region = regions{r};
    npCellsSelectivity = []; inds =[];
    i_npCellSelectivity = [];
    
    load([topDir,'AnalysesAcrossAnimals\npCells_air_',region,'.mat']);
    
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        %Remove existing tags
        filt = ['(strcmp($area, ''',region,''')) && (~isempty($SI) || (strcmp($selectivity, ''leftSelective'')) || (strcmp($selectivity, ''rightSelective'')))'];
        old = evaluatefilter(cellinfo,filt);
        
        for f = 1:size(old,1)
            cell = old(f,:);
            if isfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)},'selectivity')
                cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)} = rmfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)}, 'selectivity');
            end
            if isfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)},'SI')
                cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)} = rmfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)}, 'SI');
            end
        end
        
        %get odor-responsive cells for the animal
        cells = npCells(npCells(:,1) ==  a,[2,3,4]);

        days = unique(cells(:,1));
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            try
            load([animDir,animal,'spikes',daystr,'.mat'])
            catch
                continue
            end
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            load([animDir, animal,'nosepokeWindow',daystr,'.mat'])
            runeps = cs_getRunEpochs(animDir, animal, 'odorplace',day);
            runeps = runeps(:,2);
            
            
            %find odor triggers for the whole day
            trigs_all= []; correct_left_all = []; correct_right_all = [];
            for ep = 1:length(runeps)
                epoch = runeps(ep);
                [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                %trigs = odorTriggers{day}{epoch}.allTriggers;
                npwins = nosepokeWindow{day}{epoch};
                
                correct_left_all = [correct_left_all; npwins(correct_left,:)];
                correct_right_all = [correct_right_all; npwins(correct_right,:)];
                trigs_all = [trigs_all ; npwins];
                
            end
            
            %loop over cells and get FR for right vs left odor
            frallcellsL = []; frallcellsR = [];
            for c = 1:size(daycells,1)
                
               
                leftspikes = [];
                rightspikes = [];
                leftfr = [];
                rightfr = [];
                
                cell = daycells(c,:);
          
                 if a ==3 && day == 3 && cell(1) == 30 && cell(2) == 1
                     disp('found cell')
                 end
                runspikes = [];
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                        runspikes = [runspikes; spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1)];
                    end
                    
                    
                end
                
                for t = 1:length(correct_left_all)
                    trigwin = [correct_left_all(t,1), correct_left_all(t,2)];
                    winspikes =sum(runspikes > trigwin(1) & runspikes <= trigwin(2));
                    
                    fr = winspikes/(trigwin(2) - trigwin(1));
                    leftspikes = [leftspikes; winspikes];
                    leftfr = [leftfr; fr];
                end
                
                for t = 1:length(correct_right_all)
                    trigwin = [correct_right_all(t,1), correct_right_all(t,2)];
                    winspikes = sum(runspikes > trigwin(1) & runspikes <= trigwin(2));
                    fr = winspikes/(trigwin(2) - trigwin(1));
                    rightspikes = [rightspikes; winspikes];
                    rightfr = [rightfr; fr];
                end
                
                
                sumspikes = sum([leftspikes; rightspikes]);
                
                %number of spikes should be greater than number of trials
                %(this filter was done for finding NP cells- should no
                %exclude anything)
                if sumspikes >= length(trigs_all)
                    meanleftfr = mean(leftfr);
                    meanrightfr = mean(rightfr);
                    
                    
                   %get selectivity index
                    SI = (meanleftfr - meanrightfr)/(meanleftfr + meanrightfr);

                    %do the shuffle
                    
                    leftones = ones(size(leftspikes,1),1);
                    rightzeros = zeros(size(rightspikes,1),1);
                    
                    indicators = [leftones;rightzeros];
                    allTrials = [leftfr; rightfr];
                    
                    shuffIndexDist = zeros(iterations,1);
                    for i = 1:iterations
                        shuffledRL = indicators(randperm(length(indicators)));
                        newLeft = mean(allTrials(shuffledRL == 1));
                        newRight = mean(allTrials(shuffledRL == 0));
                        
                        newIndex = (newLeft - newRight)/(newLeft + newRight);
                        shuffIndexDist(i) = newIndex;
                    end
                    
                    shuffMean = mean(shuffIndexDist);
                    shuffStd = std(shuffIndexDist);
                    
                    %get zscore compared to shuffled distribution
                    Z = (SI - shuffMean)/shuffStd;
                    
                    
                    %cell is selective if abs val of Z > 1.5
                    for e = 1:length(runeps)
                        epoch = runeps(e);
                        if length(cellinfo{cell(1)}{epoch}{cell(2)}) >= cell(3)
                            if Z >= 1.5 
                                cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.selectivity = 'leftSelective';
                                
                            elseif Z <= -1.5
                                cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.selectivity = 'rightSelective';
                                
                            end
                            cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.SI = SI;
                        end
                    end
                end
                
            end
        end
        save([topDir,animal,'Expt\',animal,'_direct\',animal,'cellinfo.mat'],'cellinfo');
    end
end

%create matrix of cell IDs
cs_listSelectiveCells_air