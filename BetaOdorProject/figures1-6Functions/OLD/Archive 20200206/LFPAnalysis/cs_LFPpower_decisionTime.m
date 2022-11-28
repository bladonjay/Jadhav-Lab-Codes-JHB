
%Do this for each session. Get beta power on all trials, look at top and
%bottom 25% of trials.
%Get SI for these trials, for each cell
%compare SI for high and low beta trials

clear
close all
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%regions = {'CA1','PFC'};
regions = {'CA1','PFC','OB'};
fract = 0.25; %fraction of trials to take (i.e. top and bottom 25%)
freq = 'resp';

    figure, hold on
for r = 1:length(regions)
    region = regions{r};
    %for r = 1:length(regions)
        %region = regions{r};
         decisionTime_all = [];
         power_all = [];
        decisionTime_high = [];
        decisionTime_low = [];
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
                eps = cs_getRunEpochs(animDir, animal, 'odorplace',day);
                eps = eps(:,2);
                
                leftTrials = []; rightTrials = []; 
                allTrials = [];
                power = [];
                for ep = eps'
                    [cl, cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    c_left = nosepokeWindow{day}{ep}(cl,:);
                    leftTrials = [leftTrials; c_left];
                    c_right = nosepokeWindow{day}{ep}(cr,:);
                    rightTrials = [rightTrials; c_right];
                    allTrials = [allTrials; nosepokeWindow{day}{ep}(sort([cr;cl]),:)];
                    
                    %geta high vs low  trials
                    
                    %get eeg tet
                    lfptet = cs_getMostCellsTet(animal,day,ep,region);
                    
                    lfp = loadeegstruct(animDir, animal, freq,day,ep,lfptet);
                    lfp = lfp{day}{ep}{lfptet};
                    times = geteegtimes(lfp);
                    lfp = double(lfp.data(:,3));
                    
                    %get zscore
                    lfp = (lfp - mean(lfp))/std(lfp);
                    
                    lfpbins = periodAssign(times, [c_right;c_left]); %Assign spikes to align with each trials(same number = same trial, number indicates trial)
                    goodeeg = lfp(find(lfpbins));
                    lfpbins = nonzeros(lfpbins);
                    bp = [];
                    for s = unique(lfpbins)'
                        binpower = mean(goodeeg(lfpbins == s));
                        bp(s) = binpower;
                    end
                    power = [power,bp];
                end
                
                 power_all = [power_all;power'];
                decisionTime_all = [decisionTime_all; allTrials(:,2)-allTrials(:,1)];
                
                %check for np error trials
                test = find(allTrials(:,2)-allTrials(:,1) < 0.5);
                if ~isempty(test)
                    keyboard
                end
                [~,ind] = sort(power);
                allTrials = allTrials(ind,:);
               
                
                numtrials = round(fract*size(allTrials,1));   
                % --- low 
               trialTimes = allTrials(1:numtrials,:);
                exitTimes = trialTimes(:,2)-trialTimes(:,1);
                
                decisionTime_low = [decisionTime_low; exitTimes];
                
                
                % --- high
                trstart = size(allTrials,1) - numtrials +1;
               trialTimes = allTrials(trstart:end,:);
                exitTimes = trialTimes(:,2)-trialTimes(:,1);
                
                decisionTime_high = [decisionTime_high; exitTimes];
                
                
            end
            
        end

        %remove highest point from both arrays (outliers)
        [~,ind] = max(power_all);
        power_all(ind) = [];
        decisionTime_all(ind) = [];
        
        [~,ind] = max(decisionTime_all);
        power_all(ind) = [];
        decisionTime_all(ind) = [];
        
subplot(1,3,r)
plot(power_all, decisionTime_all, 'k.','MarkerSize',12)
box off
hold on

%fit = polyfit(allDivTimes(keepinds),allNPoff(keepinds),1);
fit = polyfit(power_all, decisionTime_all,1);
plot([min(power_all), max(power_all)], polyval(fit,[min(power_all), max(power_all)] ))
[CC,p] = corrcoef(power_all,decisionTime_all);

text(1, 2, ['p = ',num2str(round(p(1,2),2,'significant'))])

xlabel(['Zscored ',freq,' power'])
ylabel(['Decision Time (seconds)'])
title(region);
end

figfile = [figDir, 'EEG\',freq,'Power-DecisionTime'];
    print('-djpeg',figfile)
    print('-dpdf',figfile)