clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
lfpregions = {'CA1-PFC','CA1-OB','PFC-OB'};
freqs = {'beta','resp'};
for d = 1:length(freqs)
    
freq = freqs{d};
switch freq
    case 'beta'
        bandpass = [15 30];
    case 'resp'
        bandpass = [7 8];
end

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
                daystr = getTwoDigitNumber(day);
                epochs = dayeps(dayeps(:,1)==day,2);
                load([animDir, animal,'coherence',lfpregion,daystr]);
                
                for ep = epochs'
                    
                    disp(['Doing ',animal,' day ',num2str(day),' epoch ', num2str(ep)])
                    
                    times = coherence{day}{ep}.time;
                    data = coherence{day}{ep}.Coh;
                    goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
                    data = data(find(goodrows),:);
                    data = mean(data,1);
                    
                    
                    
                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    correctinds = [correct_left;correct_right];
                    
                    
                    windows = nosepokeWindow{day}{ep}(correctinds,:);
                    npoff = windows(:,2)-windows(:,1);
                    
                    
                    windows = [windows(:,1)-0.25, windows(:,2)+0.25];
                    peaktimes = [];
                    for w = 1:size(windows,1)
                        trialdata = data(isExcluded(times,windows(w,:)));
                        
                        if isempty(trialdata)
                            peaktime = NaN;
                        else
                        trialtimes = times(isExcluded(times,windows(w,:)));
                                   
                                    %smooth
                                    binsize = (windows(w,2)-windows(w,1))/length(trialdata);
                                    newbins = windows(w,1):binsize/3:windows(w,2)-binsize/3;
                                    trialdata = smoothdata(interp1(trialtimes,trialdata,newbins),'gaussian',15)';
                                    
                                    
                                   
                                    %look for local maxima
                                    locmax = find(islocalmax(trialdata));
                                    first = trialdata(1);
                                    last = trialdata(end);
                                    if max(trialdata(locmax)) < first || max(trialdata(locmax)) < last
                                        [~,ind] = max(trialdata);
                                    else 
                                        [~,i] = max(trialdata(locmax));
                                        ind = locmax(i);
                                    end 
                                    peaktime = newbins(ind)-windows(w,1);
                    
                        end
                        peaktimes = [peaktimes;peaktime];
                        
                        
                    end
                   
                    
                    if length(npoff) ~= length(peaktimes)
                        keyboard
                    end
%                     %remove any trials where np wasn't long enough
%                     test = npoff<0.5;
%                     if sum(test) >0
%                        npoff = npoff(~test);
%                        divTimes = divTimes(~test);
%                     end
                    
                    
                    allPeakTimes = [allPeakTimes; mean(peaktimes)];
                    allNPTimes = [allNPTimes; mean(npoff)];
                    %allPeakTimes = [allPeakTimes; peakTimes];
                end
            end
        
    end

% 
keepinds = ~isnan(allPeakTimes);
allNPTimes = allNPTimes(keepinds);
allPeakTimes = allPeakTimes(keepinds);

subplot(1,3,lr) 

plot(allNPTimes, allPeakTimes, 'k.','MarkerSize',12)
box off
hold on

fit = polyfit(allNPTimes,allPeakTimes,1);
plot([min(allNPTimes), max(allNPTimes)], polyval(fit,[min(allNPTimes), max(allNPTimes)] ))

%set(gcf,'Position',[1200,300,510,420]);
xlabel('Decision time')
ylabel(['Peak ',freq,' coherence time'])
title([lfpregion,' peak ',freq] )
[CCdiv,p_div] = corrcoef(allNPTimes,allPeakTimes);
R = CCdiv(1,2);
p = round(p_div(1,2),2,'significant');

text(1,2,['p = ', num2str(p)])

drawnow

% f = figure;
% histogram(allPeakTimes)
% pause
% close(f)
    end
    figtitle = ['DecisionTime',freq,'CoherencePeak'];
figfile = [figDir,'EEG\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);


end