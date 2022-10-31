clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
lfpregions = {'CA1-PFC','CA1-OB','PFC-OB'};

figure
    for lr = 1:length(lfpregions)
        
        
        lfpregion = lfpregions{lr};
        disp(['Starting ',lfpregion])
    
    allPeakTimes_RR = [];
    allPeakTimes_beta = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        
            nosepokeWindow= loaddatastruct(animDir,animal,'nosepokeWindow');
            odorTriggers= loaddatastruct(animDir,animal,'odorTriggers');
            tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
            
            dayeps = cs_getRunEpochs(animDir, animal,'odorplace');
            days = unique(dayeps(:,1));
            for day = days'
                daystr = getTwoDigitNumber(day);
                epochs = dayeps(dayeps(:,1)==day,2);
                load([animDir, animal,'coherence',lfpregion,daystr]);
                
                for ep = epochs'
                    
                    disp(['Doing ',animal,' day ',num2str(day),' epoch ', num2str(ep)])
                    
                    %% beta
                    bandpass = [15 30];

                    times = coherence{day}{ep}.time;
                    data = coherence{day}{ep}.Coh;
                    goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
                    data = data(find(goodrows),:);
                    data = mean(data,1);

                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    correctinds = [correct_left;correct_right];
                    
                    
                    windows = nosepokeWindow{day}{ep}(correctinds,:);
                    buffwindows = [windows(:,1) windows(:,1)+2];
                    peaktimes = [];
                    for w = 1:size(buffwindows,1)
                        triallfp = data(isExcluded(times,buffwindows(w,:)));
                        
                        if isempty(triallfp)
                            peaktime = NaN;
                        else
                        trialtimes = times(isExcluded(times,buffwindows(w,:)));
                        locmin = find(islocalmin(triallfp));
                        
                        %take first peak after the first local minimum
                        %(i.e. sometimes rhythm is high just before np so the "peak" will be the first bin, but
                        %this is not representative of the peak that occurs
                        %during sniffing). 
                        if length(locmin)>1
                        [~,ind] = max(triallfp(locmin(1):locmin(end)));
                        ind = ind + locmin(1)-1;

                        else
                            try
                            [~,ind] = max(triallfp(locmin(1):end));
                            ind = ind + locmin(1)-1;

                            catch
                            [~,ind] = max(triallfp);
                            end
                        end
%                        
                        peaktime = trialtimes(ind)-windows(w,1);
                        end
                        peaktimes = [peaktimes;peaktime];
                        
                        
                    end

                    allPeakTimes_beta = [allPeakTimes_beta; peaktimes];
                   
                    
                    %% RR 
                    bandpass = [7 8];
                    times = coherence{day}{ep}.time;
                    data = coherence{day}{ep}.Coh;
                    goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
                    data = data(find(goodrows),:);
                    data = mean(data,1);

                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    correctinds = [correct_left;correct_right];
                    
                    
                    windows = nosepokeWindow{day}{ep}(correctinds,:);
                    buffwindows = [windows(:,1) windows(:,1)+2];
                    peaktimes = [];
                    for w = 1:size(buffwindows,1)
                        triallfp = data(isExcluded(times,buffwindows(w,:)));
                        
                        if isempty(triallfp)
                            peaktime = NaN;
                        else
                        trialtimes = times(isExcluded(times,buffwindows(w,:)));
                        locmin = find(islocalmin(triallfp));
                        
                        %take first peak after the first local minimum
                        %(i.e. sometimes rhythm is high just before np so the "peak" will be the first bin, but
                        %this is not representative of the peak that occurs
                        %during sniffing). 
                        if length(locmin)>1
                        [~,ind] = max(triallfp(locmin(1):locmin(end)));
                        ind = ind + locmin(1)-1;

                        else
                            try
                            [~,ind] = max(triallfp(locmin(1):end));
                            ind = ind + locmin(1)-1;

                            catch
                            [~,ind] = max(triallfp);
                            end
                        end
%                        
                        peaktime = trialtimes(ind)-windows(w,1);
                        end
                        peaktimes = [peaktimes;peaktime];
                        
                        
                    end

                    allPeakTimes_RR = [allPeakTimes_RR; peaktimes];
                end
            end
        
    end

% 
keepinds = ~isnan(allPeakTimes_beta);
allPeakTimes_beta = allPeakTimes_beta(keepinds);
allPeakTimes_RR = allPeakTimes_RR(keepinds);

subplot(1,3,lr) 
errorbar([1 2],[mean(allPeakTimes_beta), mean(allPeakTimes_RR)], [stderr(allPeakTimes_beta), stderr(allPeakTimes_RR)],'k.');
p = ranksum(allPeakTimes_beta,allPeakTimes_RR);
box off
hold on

xlim([0 3])
ylim([0 1.5])
xticks([1 2])
xticklabels({'beta','RR'})
ylabel('Peak Coherence Time');

text(1,0.25,['p = ', num2str(p)])

drawnow

% f = figure;
% histogram(allPeakTimes)
% pause
% close(f)
    end
    figtitle = 'CoherencePeakTime';
figfile = [figDir,'EEG\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);