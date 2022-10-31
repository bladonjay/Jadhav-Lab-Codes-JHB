
clear
close all
animal = 'CS44';
day = 6; %day when learning crossed threshold
trialtype = 'novelodor';
region = 'PFC';
window = 50; %number of trials
stepsize = 1;

[topDir, figDir] = cs_setPaths;
animDir = [topDir, animal, 'Expt\',animal,'_direct\'];

% concatenate all trial times and r/l assignments for the day
odorTriggers = loaddatastruct(animDir, animal,'odorTriggers',day);
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

cellfilt = 'strcmp($selectivity, ''novelLeftSelective_postlearn'') || strcmp($selectivity, ''novelRightSelective_postlearn'')';
cells = evaluatefilter(cellinfo{day},cellfilt);

cells = unique(cells(:,[2 3]),'rows');

for c = 4:size(cells,1)
    cell = cells(c,:);
    
    goodeps = cs_findGoodEpochs(spikes{day},{'data'},cell);
    eps = intersect(goodeps, runeps);
    
    cateps = cat( 1, spikes{day}(eps'));
    allSpikes = cellfun(@(x) x{cell(1)}{cell(2)}.data(:,1), cateps,'UniformOutput',false);
    allSpikes = cat(1,allSpikes{:});
    
    t = 1;
    selectivity_allwins = [];
    %get spiking within np for set number of trials
    while t + window -1 < size(allTrials,1)
        wins = allTrials(t:t+window,:);
        
        leftfr = []; rightfr = [];
        for w = 1:size(wins,1)
            winspikes = isExcluded(allSpikes,wins(w,:));
            fr =sum(winspikes)/(wins(w,2) - wins(w,1));
            
            if ismember(wins(w,:), leftTrials,'rows')
                leftfr = [leftfr; fr];
            elseif ismember(wins(w,:), rightTrials,'rows')
                rightfr = [rightfr; fr];
            else
                warning('Cannot match trial window to right or left')
            end
        end
        %calculate selectivity index
        meanleftfr = mean(leftfr);
        meanrightfr = mean(rightfr);
        selectivity = (meanleftfr - meanrightfr)/(meanleftfr + meanrightfr);
        selectivity_allwins = [selectivity_allwins; selectivity];
        
        %move window and repeat
        t = t + stepsize;
    end
    
    selectivity_allwins = abs(selectivity_allwins);
    figure,
    plot(selectivity_allwins);
    
    
    
  
    
    xlabel('Trial Window')
    ylabel('Selectivity Index Magnitude');
    axis([0 length(selectivity_allwins) 0 max(selectivity_allwins) + 0.1*max(selectivity_allwins)])
    
    figfile = [figDir,'OdorSelectivity\',animal,'_',trialtype,'_',num2str(day),'_',num2str(cell(1)),'_',num2str(cell(2))];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        close all
    
end



