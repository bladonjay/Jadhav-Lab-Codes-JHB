%function cs_selectivityTag(cellinfo, cellind, spikes, lefttrigs, righttrigs, win, iterations)

%combine spiketimes across run epochs
%combine trig times across run epochs

%tag cells with selectivity index #, then can find these cells again and
%plot.
topDir = cs_setPaths;

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS39'};
regions = {'CA1','PFC'};
%win = [0 1];
iterations = 1000;

%winsize = win(2) + win(1);

for r = 1:length(regions)
    region = regions{r};
    npCellsSelectivity = []; inds =[];
    i_npCellSelectivity = [];
    total = 0;
    
    load([topDir,'AnalysesAcrossAnimals\npCells_',region,'.mat']);
    
    for a = 8 %1:length(animals)
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
        
        %         filt = ['(~isempty($SI))'];
        %         test = evaluatefilter(cellinfo,filt);
        %          for f = 1:size(test,1)
        %             cell = test(f,:);
        %             if isfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)},'SI')
        %                 cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)} = rmfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)}, 'SI');
        %             end
        %         end
        %
        
        % Im guessing npcells is loaded in
        cells = npCells(npCells(:,1) ==  a,[2,3,4]);
        
        days = unique(cells(:,1));
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            % daystr=sprintf('%0.2d',day); % this works better
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat']) % this is the cell spikes
            load([animDir,animal,'odorTriggers',daystr,'.mat']) % this is which odor goes to which side
            load([animDir, animal,'nosepokeWindow',daystr,'.mat']) % this is the start/stop times for the pokes
            runeps = find(~cellfun(@isempty,nosepokeWindow{day}));
            
            
            trigs_all= []; correct_left_all = []; correct_right_all = [];
            
            % this concatenates all the odor triggers
            for ep = 1:length(runeps)
                epoch = runeps(ep);
                [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                %trigs = odorTriggers{day}{epoch}.allTriggers;
                npwins = nosepokeWindow{day}{epoch};
                
                correct_left_all = [correct_left_all; npwins(correct_left,:)];
                correct_right_all = [correct_right_all; npwins(correct_right,:)];
                trigs_all = [trigs_all ; npwins];
                
            end
            
            frallcellsL = []; frallcellsR = [];
            for c = 1:size(daycells,1)
                
                leftspikes = [];
                rightspikes = [];
                leftfr = [];
                rightfr = [];
                
                cell = daycells(c,:);
          
                runspikes = [];
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        runspikes = [runspikes; spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1)];
                    end
                    
                    
                end
                
                for t = 1:length(correct_left_all)
                    %trig = correct_left_all(t);
                    trigwin = [correct_left_all(t,1), correct_left_all(t,2)];
                    % sum all spikes after the poke start and before the
                    % poke end
                    winspikes =sum(runspikes > trigwin(1) & runspikes <= trigwin(2));
                    % divide them by the poke duration
                    %bins = (trig-win(1):binsize:trig+win(2));
                    %binspikecount = histcounts(winspikes,bins);
                    fr = winspikes/(trigwin(2) - trigwin(1));
                    leftspikes = [leftspikes; winspikes];
                    leftfr = [leftfr; fr];
                end
                
                for t = 1:length(correct_right_all)
                    %trig = correct_right_all(t);
                    trigwin = [correct_right_all(t,1), correct_right_all(t,2)];
                    winspikes = sum(runspikes > trigwin(1) & runspikes <= trigwin(2));
                    %bins = (trig-win(1):binsize:trig+win(2));
                    %binspikecount = histcounts(winspikes,bins);
                    fr = winspikes/(trigwin(2) - trigwin(1));
                    rightspikes = [rightspikes; winspikes];
                    rightfr = [rightfr; fr];
                end
                
                
                sumspikes = sum([leftspikes; rightspikes]);
                
                % if there are more spikes than trials
                if sumspikes >= length(trigs_all)
                    meanleftfr = mean(leftfr);
                    meanrightfr = mean(rightfr);
                    
                    
                    selectivity = (meanleftfr - meanrightfr)/(meanleftfr + meanrightfr);
                    
                    
                    
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
                    % so this is the zscore of the null
                    trueIndexZ = (selectivity - shuffMean)/shuffStd;
                    
                    % for each epoch, does the z exceed 1.5 sd over the
                    % mean... so a p of 0.0668
                    for e = 1:length(runeps)
                        epoch = runeps(e);
                        if length(cellinfo{cell(1)}{epoch}{cell(2)}) >= cell(3)
                            if trueIndexZ >= 1.5
                                cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.selectivity = 'leftSelective';
                                total = total+1;
                            elseif trueIndexZ <= -1.5
                                cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.selectivity = 'rightSelective';
                                total = total+1;
                            end
                            cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.SI = selectivity;
                        end
                    end
                end
                
            end
        end
        save([topDir,animal,'Expt\',animal,'_direct\',animal,'cellinfo.mat'],'cellinfo');
    end
end

cs_listSelectiveCells