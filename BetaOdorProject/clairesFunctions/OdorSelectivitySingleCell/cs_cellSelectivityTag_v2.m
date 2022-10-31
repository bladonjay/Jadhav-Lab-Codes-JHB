%only look at spikes within high beta times

clear
topDir = cs_setPaths;

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS34'};
regions = {'CA1','PFC'};
%win = [0 1];
iterations = 1000;

%winsize = win(2) + win(1);

for r = 1:length(regions)
    selectiveCT=0;
    region = regions{r};
    npCellsSelectivity = []; inds =[];
    i_npCellSelectivity = [];
    load([topDir,'AnalysesAcrossAnimals\npCells_',region,'.mat']);
    
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        load([animDir, animal,'betaWindows.mat']) %this should give start/end time for beta on each trial
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
        cells = npCells(npCells(:,1) ==  a,[2,3,4]);
        
        %         filt = ['isequal($area,''',region,''') && strcmp($tag, ''accepted'')'];
        %         singleunit = evaluatefilter(cellinfo,filt);
        %
        %         cells = intersect(cells, singleunit(:,[1 3 4]),'rows');
        
        days = unique(cells(:,1));
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            %load([animDir, animal,'nosepokeWindow',daystr,'.mat'])
            %runeps = find(~cellfun(@isempty,nosepokeWindow{day}));
            runeps = cs_getRunEpochs(animDir,animal,'odorplace',day);
            runeps = runeps(:,2);
            
            trigs_all= []; correct_left_all = []; correct_right_all = [];
            for ep = 1:length(runeps)
                epoch = runeps(ep);
                [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                %trigs = odorTriggers{day}{epoch}.allTriggers;
                betaWinsAll = betaWindows{day}{epoch};
                betaWins = betaWinsAll(correct_left,:);
                betaWins(~any(~isnan(betaWins), 2),:)=[];
                correct_left_all = [correct_left_all; betaWins];
                
                betaWins = betaWinsAll(correct_right,:);
                betaWins(~any(~isnan(betaWins), 2),:)=[];
                correct_right_all = [correct_right_all; betaWins];
                trigs_all = [trigs_all; betaWinsAll];
                
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
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
%                         if strcmp(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.tag, 'mua')
%                             continue
%                         else
                            runspikes = [runspikes; spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1)];
%                         end
                    end
                    
                    
                end
                
                for t = 1:size(correct_left_all,1)
                    %trig = correct_left_all(t);
                    trigwin = [correct_left_all(t,1), correct_left_all(t,2)];
                    winspikes =sum(runspikes > trigwin(1) & runspikes <= trigwin(2));
                    %bins = (trig-win(1):binsize:trig+win(2));
                    %binspikecount = histcounts(winspikes,bins);
                    fr = winspikes/(trigwin(2) - trigwin(1));
                    leftspikes = [leftspikes; winspikes];
                    leftfr = [leftfr; fr];
                end
                
                for t = 1:size(correct_right_all,1)
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
                    
                    trueIndexZ = (selectivity - shuffMean)/shuffStd;
                    
                    
                    
                    for e = 1:length(runeps)
                        epoch = runeps(e);
                        if length(cellinfo{cell(1)}{epoch}{cell(2)}) >= cell(3)
                            if trueIndexZ >= 1.5
                                cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.selectivity = 'leftSelective';
                                selectiveCT=selectiveCT+1;
                            elseif trueIndexZ <= -1.5
                                cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.selectivity = 'rightSelective';
                                selectiveCT=selectiveCT+1;
                            end
                            cellinfo{cell(1)}{epoch}{cell(2)}{cell(3)}.SI = selectivity;
                        end
                    end
                    
                end
            end
        end
        save([topDir,animal,'Expt\',animal,'_direct\',animal,'cellinfo.mat'],'cellinfo');
    end
    fprintf('in %s there are %d selective cells \n',regions{r},selectiveCT);
end

cs_listSelectiveCells