
clear
close all
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%regions = {'CA1','PFC'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};
fract = 0.25;
freq = 'resp';

switch freq
    case 'beta'
        bandpass = [15 30];
    case 'resp'
        bandpass = [7 8];
end

    figure, hold on
for r = 1:length(regions)
    region = regions{r};
    %for r = 1:length(regions)
        %region = regions{r};
        decisionTime_all = [];
         coh_all = [];
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
            %tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
            
            for day = days'
                daystr = getTwoDigitNumber(day);
                nosepokeWindow = loaddatastruct(animDir, animal,'nosepokeWindow',day);
                eps = cs_getRunEpochs(animDir, animal, 'odorplace',day);
                eps = eps(:,2);
                
                leftTrials = []; rightTrials = []; allTrials = [];
                
                load([animDir, animal,'coherence',region,daystr]);
                
                coh = [];
                for ep = eps'
                    [cl, cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    c_left = nosepokeWindow{day}{ep}(cl,:);
                    leftTrials = [leftTrials; c_left];
                    c_right = nosepokeWindow{day}{ep}(cr,:);
                    rightTrials = [rightTrials; c_right];
                    allTrials = [allTrials; nosepokeWindow{day}{ep}(sort([cr;cl]),:)];
                    
                    %geta high vs low  trials

                    
                    times = coherence{day}{ep}.time;
                    data = coherence{day}{ep}.Coh;
                    goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
                    data = data(find(goodrows),:);
                    data = mean(data,1);
                    
                    trialtimes = sort([c_right;c_left]);
                    lfpbins = periodAssign(times, trialtimes); %Assign spikes to align with each trials(same number = same trial, number indicates trial)
                      if length(unique(nonzeros(lfpbins))) ~= size([c_right;c_left],1)       
                          allTrials(end,:)  = [];
                      end
                      
                    goodeeg = data(find(lfpbins));
                    lfpbins = nonzeros(lfpbins);
                    bp = [];
                    for s = unique(lfpbins)'
                        binpower = mean(goodeeg(lfpbins == s));
                       
                        bp(s) = binpower;
                    end
                    coh = [coh,bp];
                end
                
                coh_all = [coh_all;coh'];
                decisionTime_all = [decisionTime_all; allTrials(:,2)-allTrials(:,1)];
                
               
                [coh,ind] = sort(coh);
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
         
        
subplot(3,1,r)
plot(coh_all, decisionTime_all, 'k.','MarkerSize',12)
box off
hold on

%fit = polyfit(allDivTimes(keepinds),allNPoff(keepinds),1);
fit = polyfit(coh_all, decisionTime_all,1);
plot([min(coh_all), max(coh_all)], polyval(fit,[min(coh_all), max(coh_all)] ))
[CC,p] = corrcoef(coh_all,decisionTime_all);

text(1, 2, ['p = ',num2str(round(p(1,2),2,'significant'))])

xlabel(['Zscored ',freq,' coherence'])
ylabel(['Decision Time (seconds)'])
title(region);
    

%     mn = [mean(decisionTime_low),mean(decisionTime_high)];
%     stdev = [std(decisionTime_low), std(decisionTime_high)];
%     err = [stderr(decisionTime_low),stderr(decisionTime_high)];
%     
%     errorbar(r,mn(1),err(1),'r.')
%     errorbar(r+0.5,mn(2),err(2),'b.')
% 
%     %axis([0 3 0 (max(mn)+2*max(stdev))])
%     [p] = ranksum(decisionTime_low,decisionTime_high);
%     text(r, (max(mn)+max(err)), ['p = ',num2str(round(p,2))])
    
end
%axis([0 4 0 1]);
%xlim([0 4])
% xticks(1:length(regions))
%     xticklabels(regions);
%     ylabel('Decision Time')
    
    figfile = [figDir, 'EEG\',freq,'Coherence-DecisionTime'];
    print('-djpeg',figfile)
    print('-dpdf',figfile)