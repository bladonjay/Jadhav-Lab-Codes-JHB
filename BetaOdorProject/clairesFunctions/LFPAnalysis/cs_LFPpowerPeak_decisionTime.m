clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
lfpregions = {'CA1','PFC','OB'};
freqs = {'beta','resp'};
for d = 1:length(freqs)
    
freq = freqs{d};

figure
    for lr = 1:length(lfpregions)
        
        
        lfpregion = lfpregions{lr};
        disp(['Starting ',lfpregion])
    allNPTimes = [];
    allPeakTimes = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        
            nosepokeWindow= loaddatastruct(animDir,animal,'nosepokeWindow');
            odorTriggers= loaddatastruct(animDir,animal,'odorTriggers');
            tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
            
            dayeps = cs_getRunEpochs(animDir, animal,'odorplace');
            days = unique(dayeps(:,1));
            for day = days'
                epochs = dayeps(dayeps(:,1)==day,2);
                for ep = epochs'
                    
                    disp(['Doing ',animal,' day ',num2str(day),' epoch ', num2str(ep)])
                    lfptet = cs_getMostCellsTet(animal,day,ep,lfpregion);
                    
                    lfp = loadeegstruct(animDir, animal, freq,day,ep,lfptet);
                    lfp = lfp{day}{ep}{lfptet};
                    times = geteegtimes(lfp);
                    lfp = double(lfp.data(:,3));
                    
                    
                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    correctinds = [correct_left;correct_right];
                    
                    
                    windows = nosepokeWindow{day}{ep}(correctinds,:);
                    buffwindows = [windows(:,1) windows(:,1)+2];
                    peaktimes = [];
                    for w = 1:size(buffwindows,1)
                        triallfp = lfp(isExcluded(times,buffwindows(w,:)));
                        
                        if isempty(triallfp)
                            peaktime = NaN;
                        else
                        trialtimes = times(isExcluded(times,buffwindows(w,:)));
                        locmin = find(islocalmin(triallfp));
                        
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
                    npoff = windows(:,2)-windows(:,1);
                    
                    if length(npoff) ~= length(peaktimes)
                        keyboard
                    end
%                     %remove any trials where np wasn't long enough
%                     test = npoff<0.5;
%                     if sum(test) >0
%                        npoff = npoff(~test);
%                        divTimes = divTimes(~test);
%                     end
                    
                    
                    allPeakTimes = [allPeakTimes; peaktimes];
                    allNPTimes = [allNPTimes; npoff];
                    %allPeakTimes = [allPeakTimes; peakTimes];
                end
            end
        
    end

% 
keepinds = ~isnan(allPeakTimes);
allNPTimes = allNPTimes(keepinds);
allPeakTimes = allPeakTimes(keepinds);

subplot(1,3,lr) 
drawnow
plot(allNPTimes, allPeakTimes, 'k.','MarkerSize',12)
box off
hold on

fit = polyfit(allNPTimes,allPeakTimes,1);
plot([min(allNPTimes), max(allNPTimes)], polyval(fit,[min(allNPTimes), max(allNPTimes)] ))

set(gcf,'Position',[1200,300,510,420]);
xlabel('Decision time')
ylabel(['Peak ',freq,' power time'])
title([lfpregion,' peak ',freq] )
[CCdiv,p_div] = corrcoef(allNPTimes,allPeakTimes);
R = CCdiv(1,2);
p = p_div(1,2);

text(1,1,['R = ',num2str(R), newline, 'p = ', num2str(p)])



    end
    figtitle = ['DecisionTime',freq,'PowerPeak'];
figfile = [figDir,'EEG\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);
end