
clear
close all
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%regions = {'CA1','PFC'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};

freqs = {'beta','resp'};
fignum = 0;
figure, hold on
for f = 1:length(freqs)
    freq = freqs{f};
    
    switch freq
        case 'beta'
            bandpass = [15 30];
        case 'resp'
            bandpass = [7 8];
    end
    
    
    for r = 1:length(regions)
        region = regions{r};
        fignum = fignum+1;
        %for r = 1:length(regions)
        %region = regions{r};
        decisionTime_all = [];
        coh_all = [];
        
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            disp(['Doing ',animal,' ',freq,' ',region]);
            runeps = cs_getRunEpochs(animDir, animal, 'odorplace');
            if isempty(runeps)
                continue
            end
            days = unique(runeps(:,1));
            rewards = loaddatastruct(animDir, animal,'rewards');
            odorTriggers = loaddatastruct(animDir, animal,'odorTriggers',days);
            tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
            
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
                    switch freq
                        case 'beta'
                            
                            data = coherence{day}{ep}.rawCoh;
                            sd = coherence{day}{ep}.rawSD;
                            mn = coherence{day}{ep}.rawMean;
                            data = (data-mn)./sd;
                        case 'resp'
                            
                            data = coherence{day}{ep}.rawCoh;
                            time = coherence{day}{ep}.time;
                            %get baseline for zscore
                            rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
                            rewarddata = data(:,isExcluded(time,rtimes));
                            mn = mean(rewarddata,2);
                            sd = std(rewarddata,0,2);
                
                            %zscore to reward coh
                            data = (data-mn)./sd;
                    end
                    
                    goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
                    data = data(find(goodrows),:);
                    data = mean(data,1);
                    
                    trialtimes = sort([c_right;c_left]);
                    lfpbins = periodAssign(times, trialtimes); %Assign data to align with each trials(same number = same trial, number indicates trial)
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
                
                
%                 [coh,ind] = sort(coh);
%                 allTrials = allTrials(ind,:);
%                 
%                 numtrials = round(fract*size(allTrials,1));
                
                
                
%                 % --- low
%                 trialTimes = allTrials(1:numtrials,:);
%                 exitTimes = trialTimes(:,2)-trialTimes(:,1);
%                 
%                 decisionTime_low = [decisionTime_low; exitTimes];
%                 
%                 
%                 % --- high
%                 trstart = size(allTrials,1) - numtrials +1;
%                 trialTimes = allTrials(trstart:end,:);
%                 exitTimes = trialTimes(:,2)-trialTimes(:,1);
%                 
%                 decisionTime_high = [decisionTime_high; exitTimes];
                
                
            end
            
        end
       
        
        [CC,p] = corrcoef(coh_all,decisionTime_all);
        p = p(1,2);
        R = CC(1,2);
        %get sextiles
        
        subplot(2,3,fignum)
        
        [means, errs] = cs_sextiles(decisionTime_all, coh_all);

       
        ylabel([freq, ' coherence'])
        xlabel(['Decision Latency'])
        %text(1,means(1),['p = ',num2str(round(p,2,'significant'))])
        %plot(coh_all, decisionTime_all, 'k.','MarkerSize',12)
        box off
        hold on
        text(1, means(1), ['R = ',num2str(round(R,2,'significant')),newline,'p = ',num2str(round(p,2,'significant'))  ])
        title(region)
        drawnow
        
        
    end
    
end

figfile = [figDir, 'EEG\Coherence-DecisionTime'];
print('-djpeg',figfile)
print('-dpdf',figfile)