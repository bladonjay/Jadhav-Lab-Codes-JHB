%SI index on high vs low beta trials

%Do this for each session. Get beta power on all trials, look at top and
%bottom 25% of trials.
%Get SI for these trials, for each cell
%compare SI for high and low beta trials

clear
close all
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
betaregions = {'CA1','PFC','OB'};
numtrials = 10;
freq = 'resp';

[A,B] = meshgrid(1:length(cellregions),1:length(betaregions));
pairInds=reshape(cat(2,A',B'),[],2);

figure
hold on
for r = 1:size(pairInds,1)
    cellregion = cellregions{pairInds(r,1)};
    betaregion = betaregions{pairInds(r,2)};
    
    selectivity_high = [];
    selectivity_low = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        runeps = cs_getRunEpochs(animDir, animal, 'odorplace');
        if isempty(runeps)
            continue
        end
        days = unique(runeps(:,1));
        
        odorTriggers = loaddatastruct(animDir, animal,'odorTriggers',days);
        tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
        
        for day = days'
            nosepokeWindow = loaddatastruct(animDir, animal,'nosepokeWindow',day);
            spikes = loaddatastruct(animDir, animal,'spikes',day);
            eps = cs_getRunEpochs(animDir, animal, 'odorplace',day);
            eps = eps(:,2);
            
            leftTrials = []; rightTrials = []; allTrials = [];
            allbeta =[];
            
            
            betapower = [];
            for ep = eps'
                [cl, cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                c_left = nosepokeWindow{day}{ep}(cl,:);
                leftTrials = [leftTrials; c_left];
                c_right = nosepokeWindow{day}{ep}(cr,:);
                rightTrials = [rightTrials; c_right];
                allTrials = [allTrials; nosepokeWindow{day}{ep}(sort([cr;cl]),:)];
                
                %geta high vs low beta trials
                
                %get beta tet
                tetfilt = ['strcmp($area,''',betaregion,''')'];
                betatet = evaluatefilter(tetinfo{day}{ep},tetfilt);
                if length(betatet) > 1
                    numcells = cellfetch(tetinfo{day}{ep}(betatet),'numcells');
                    [~, ind] = max(cell2mat(numcells.values));
                    if length(ind) > 1 %if more than one tet have highest num cells, just take the first
                        ind = ind(1);
                    end
                    betatet = betatet(ind);
                end
                
                eeg = loadeegstruct(animDir, animal, freq,day,ep,betatet);
                eeg = eeg{day}{ep}{betatet};
                times = geteegtimes(eeg);
                eeg = double(eeg.data(:,3));
                
                betabins = periodAssign(times, [c_right;c_left]); %Assign spikes to align with each trials(same number = same trial, number indicates trial)
                goodbeta = eeg(find(betabins));
                betabins = nonzeros(betabins);
                bp = [];
                for s = unique(betabins)'
                    power = mean(goodbeta(betabins == s));
                    bp(s) = power;
                end
                betapower = [betapower,bp];
            end
            
            [betapower,ind] = sort(betapower);
            allTrials = allTrials(ind,:);
            
            %get cells on that day
            cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
            
            cellfilt = ['strcmp($area,''',cellregion,''') && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective''))'];
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
                
                % --- low beta
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
                selectivity_low = [selectivity_low; selectivity];
                
                
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
                selectivity_high = [selectivity_high; selectivity];
            end
            
        end
    end
    
    
    
    mn = [mean(selectivity_low),mean(selectivity_high)];
    err = [stderr(selectivity_low),stderr(selectivity_high)];
    errorbar(r,mn(1),err(1),'r.')
    errorbar(r+0.25,mn(2),err(2),'k.')

    axis([0 7 0 0.5])
    x{r} = ([cellregion, ' cells, ',betaregion,' ',freq]);
    [h,p] = ttest(selectivity_low,selectivity_high);
    text(r, 0.05, ['p = ',num2str(round(p,2))])
    %xticks([1 2])
    
end
xticks([1:length(pairInds)])
xticklabels(x);
xtickangle(45)
ylabel('Selectivity Index Magnitude')
legend({['Low ',freq,' Power'],['High ',freq,' Power']})

figfile = [figDir, 'OdorSelectivity\Power-SelectivityIndex_',freq];
print('-djpeg',figfile)
print('-dpdf',figfile)
close