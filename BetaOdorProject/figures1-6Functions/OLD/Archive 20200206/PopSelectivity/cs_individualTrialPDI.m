clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
regions = {'CA1','PFC'};
lfpregions = {'CA1','PFC','OB'};
freq = 'resp';
for r = 1:length(regions)
    region = regions{r};
    
    figure
    for lr = 1:length(lfpregions)
        lfpregion = lfpregions{lr};
    allPeakTimes = [];
    allDivTimes = [];
    allPeakTimes = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        files = dir([animDir,animal,'trialSigPDI_',region,'*']);
        if ~isempty(files)
            
            nosepokeWindow= loaddatastruct(animDir,animal,'nosepokeWindow');
            odorTriggers= loaddatastruct(animDir,animal,'odorTriggers');
                 tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
                 
            for f = 1:length(files)
                load([animDir,files(f).name]);
                day = length(trialSigPDI);
                daystr = getTwoDigitNumber(day);
                %load([animaldir,animal,'trialSigPDI_',region,daystr,'.mat'])

                epochs = find(~cellfun(@isempty,trialSigPDI{day}));
                for e = 1:length(epochs)
                    ep = epochs(e);
                    
                    lfptet = cs_getMostCellsTet(animal,day,ep,lfpregion);
                    
                    lfp = loadeegstruct(animDir, animal, freq,day,ep,lfptet);
                    lfp = lfp{day}{ep}{lfptet};
                    times = geteegtimes(lfp);
                    lfp = double(lfp.data(:,3));
                    
                    
                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    correctinds = [correct_left;correct_right];
                    
                    
                    windows = nosepokeWindow{day}{ep}(correctinds,:);
                    peaktimes = [];
                    for w = 1:size(windows,1)
                        triallfp = lfp(isExcluded(times,windows(w,:)));
                        
                        if isempty(triallfp)
                            peaktime = NaN;
                        else
                        trialtimes = times(isExcluded(times,windows(w,:)));
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
                    divTimes = trialSigPDI{day}{ep};
                    
                    if length(divTimes) ~= length(peaktimes)
                        keyboard
                    end
%                     %remove any trials where np wasn't long enough
%                     test = npoff<0.5;
%                     if sum(test) >0
%                        npoff = npoff(~test);
%                        divTimes = divTimes(~test);
%                     end
                    
                    
                    allPeakTimes = [allPeakTimes; peaktimes];
                    allDivTimes = [allDivTimes; divTimes];
                    %allPeakTimes = [allPeakTimes; peakTimes];
                end
            end
        end
    end

% 
keepinds = ~isnan(allDivTimes(:,1));
allDivTimes = allDivTimes(keepinds);
allPeakTimes = allPeakTimes(keepinds);

subplot(1,3,lr) 
plot(allDivTimes, allPeakTimes, 'k.','MarkerSize',12)
box off
hold on

fit = polyfit(allDivTimes,allPeakTimes,1);
plot([min(allDivTimes), max(allDivTimes)], polyval(fit,[min(allDivTimes), max(allDivTimes)] ))

set(gcf,'Position',[1200,300,510,420]);
xlabel('First significant divergence time')
ylabel(['Peak ',freq,' power time'])
title([region,'PV, ',lfpregion,' peak ',freq] )
[CCdiv,p_div] = corrcoef(allDivTimes,allPeakTimes);
R = CCdiv(1,2);
p = p_div(1,2);

text(0.2,1,['R = ',num2str(R), newline, 'p = ', num2str(p)])

figtitle = ['PVDivergence',region,'-',freq,'PowerPeak',lfpregion];
figfile = [figDir,'PopSelectivity\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);

    end
end