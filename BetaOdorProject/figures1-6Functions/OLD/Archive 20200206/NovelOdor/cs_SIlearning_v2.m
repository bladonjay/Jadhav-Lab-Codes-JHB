%cs_siLearning2

%cells that eventually are selective, how does SI change over time (mean
%across cells, look at early learning vs late learning)

%look at day when learning crosses threshold for each animal


clear
close all
[topDir, figDir] = cs_setPaths;
animals = {'CS39','CS41','CS42','CS44'};
trialtypes = {'novelodor','novelodor2'};
regions = {'CA1','PFC'};
numtrials = 20; %number of trials
stepsize = 1;

for r = 1:length(regions)
    region = regions{r};
    
    selectivity_early = [];
    selectivity_late = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        for tt = 1:length(trialtypes)
            trialtype = trialtypes{tt};
            %find the day when learning took place
            
            runeps = cs_getRunEpochs(animDir, animal, trialtype);
            if isempty(runeps)
                continue
            end
            odorTriggers = loaddatastruct(animDir, animal,'odorTriggers',unique(runeps(:,1)));
            
            epfilt = '~isempty($prelearn) && ~isempty($postlearn)';
            eps = evaluatefilter(odorTriggers,epfilt);
            
            if isempty(eps)
                continue
            end
            day = eps(1);
            
            
            % concatenate all trial times and r/l assignments for the day
            
            nosepokeWindow = loaddatastruct(animDir, animal,'nosepokeWindow',day);
            spikes = loaddatastruct(animDir, animal,'spikes',day);
            runeps = cs_getRunEpochs(animDir, animal, trialtype,day);
            runeps = runeps(:,2);
            %windows = cat( 1, nosepokeWindow{day}{:});
            leftTrials = []; rightTrials = []; allTrials = [];
            for ep = runeps'
                [cl, cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                c_left = nosepokeWindow{day}{ep}(cl,:);
                leftTrials = [leftTrials; c_left];
                c_right = nosepokeWindow{day}{ep}(cr,:);
                rightTrials = [rightTrials; c_right];
                allTrials = [allTrials; nosepokeWindow{day}{ep}(sort([cr;cl]),:)];
            end
            
            
            % Find cells on day when learning crossed threshold
            cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
            
            cellfilt = ['strcmp($area,''',region,''') && (strcmp($selectivity, ''novelLeftSelective_postlearn'') || strcmp($selectivity, ''novelRightSelective_postlearn''))'];
            %cellfilt = ['strcmp($area,''',region,''')'];
            cells = evaluatefilter(cellinfo{day},cellfilt);
            
            cells = unique(cells(:,[2 3]),'rows');
            
            for c = 1:size(cells,1)
                cell = cells(c,:);
                
                goodeps = cs_findGoodEpochs(spikes{day},{'data'},cell);
                eps = intersect(goodeps, runeps);
                
                cateps = cat( 1, spikes{day}(eps'));
                allSpikes = cellfun(@(x) x{cell(1)}{cell(2)}.data(:,1), cateps,'UniformOutput',false);
                allSpikes = cat(1,allSpikes{:});
                
                % --- early learning
                wins = allTrials(1:numtrials,:);
                leftfr = []; rightfr = [];
                totspikes = 0;
                for w = 1:size(wins,1)
                    winspikes = isExcluded(allSpikes,wins(w,:));
                    
                    totspikes = totspikes+sum(winspikes);
                    
                    fr =sum(winspikes)/(wins(w,2) - wins(w,1));
                    
                    if ismember(wins(w,:), leftTrials,'rows')
                        leftfr = [leftfr; fr];
                    elseif ismember(wins(w,:), rightTrials,'rows')
                        rightfr = [rightfr; fr];
                    else
                        warning('Cannot match trial window to right or left')
                    end
                end
                meanleftfr = mean(leftfr);
                meanrightfr = mean(rightfr);
                
                %use absolute value- farther from zero = more selective.
                %does not matter which direction.
                if (meanleftfr + meanrightfr) == 0 
                    selectivity = 0;
                    
                elseif totspikes < numtrials
                    %don't use if not enough spikes
                    selectivity = 0;
                    
                else
                    selectivity = abs((meanleftfr - meanrightfr)/(meanleftfr + meanrightfr));
                end
                selectivity_early = [selectivity_early; selectivity];
                
                
                % --- late learning
                trstart = size(allTrials,1) - numtrials +1;
                wins = allTrials(trstart:end,:);
                leftfr = []; rightfr = [];
                totspikes = 0;
                for w = 1:size(wins,1)
                    winspikes = isExcluded(allSpikes,wins(w,:));
                    fr =sum(winspikes)/(wins(w,2) - wins(w,1));
                    
                    totspikes = totspikes+sum(winspikes);
                    
                    if ismember(wins(w,:), leftTrials,'rows')
                        leftfr = [leftfr; fr];
                    elseif ismember(wins(w,:), rightTrials,'rows')
                        rightfr = [rightfr; fr];
                    else
                        warning('Cannot match trial window to right or left')
                    end
                end
                meanleftfr = mean(leftfr);
                meanrightfr = mean(rightfr);
                
                if (meanleftfr + meanrightfr) == 0 
                    selectivity = 0;
                    
                elseif totspikes < numtrials
                %don't use if FR is too low
                selectivity = 0;
                
                else
                    selectivity = abs((meanleftfr - meanrightfr)/(meanleftfr + meanrightfr));
                end
                selectivity_late = [selectivity_late; selectivity];
                
            end
      
        end
    end
    
    sel = [selectivity_early, selectivity_late];
    
    %remove all rows where selectivity is zero
    bad = all(sel == 0,2);
    sel(bad,:) = [];
    
    sel = sel(all(~isnan(sel),2),:);
    
    if isempty(sel)
        continue
    end
    mnsel = mean(sel,1);
    err = stderr(sel);
    figure
    
    if size(sel,1) > 1
    errorbar(mnsel,err,'k.')
    else
        plot([1 2],sel,'k.')
    end
    xticks([1 2]);
    xticklabels({'prelearn','postlearn'});
    axis([0 3 0 max(mnsel)+(2*max(err))])
    ylabel('Selectivity Index Magnitude')
    
    p = ranksum(sel(:,1),sel(:,2));
    
    title(region)
    
    figfile = [figDir, 'OdorSelectivity\silearning_',region];
    print('-dpdf',figfile);
    print('-djpeg',figfile);
    
    
end

