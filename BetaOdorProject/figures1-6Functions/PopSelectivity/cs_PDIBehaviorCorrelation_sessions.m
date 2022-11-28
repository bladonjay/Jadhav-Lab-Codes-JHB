%raster showing timing of odor onset, neural activity divergence (PDI), and
%nosepoke exit

%fig1: raster
%fig2: correlation between PDI and NP exit
%fig3: comparison between CA1 and PFC mean digergence time, and time
%between divergence and NP exit (with STD) - may show which one is more
%reliable?

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};
%animals = {'CS33'};
regions = {'CA1','PFC'};

for r = 1:length(regions)
    region = regions{r};
    
    allNPoff = [];
    allDivTimes = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animaldir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        files = dir([animaldir,animal,'sessionSigPDI_',region,'*']);
        
        if ~isempty(files)
            for f = 1:length(files)
                load(files(f).name);
                day = length(sessionSigPDI);
                daystr = getTwoDigitNumber(day);
                load([animaldir,animal,'nosepokeWindow',daystr,'.mat'])
                load([animaldir,animal,'odorTriggers',daystr,'.mat'])
                
                %load([animaldir,animal,'trialSigPDI_',region,daystr,'.mat'])

                epochs = find(~cellfun(@isempty,sessionSigPDI{day}));
                for e = 1:length(epochs)
                    epoch = epochs(e);

                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                    correctinds = sort([correct_left;correct_right]);
                    
                    windows = nosepokeWindow{day}{epoch}(correctinds,:);
                    
                    npoff = windows(:,2)-windows(:,1); 
                    meannpoff = mean(npoff); %normalize time windows
                    divTimes = sessionSigPDI{day}{epoch};
                    %peakTimes = trialPeakPDI{day}{epoch}(correctinds,:);
                    
                    allNPoff = [allNPoff; meannpoff];
                    allDivTimes = [allDivTimes; divTimes];
                    %allPeakTimes = [allPeakTimes; peakTimes];
                end
            end
        end
    end
%% Figure 1 - raster
% [sortedPDI,ind] = sort(allDivTimes);
% sortedNP = allNPoff(ind);
% 
% figure,
% plot([zeros(length(allNPoff),1), sortedNP], [(1:length(allNPoff))', (1:length(allNPoff))'], 'rx')
% hold on
% plot(sortedPDI, 1:length(allDivTimes), 'go')



[sortedNP, ind] = sort(allNPoff);
sortedPDI = allDivTimes(ind);
%sortedPeaks = allPeakTimes(ind);
bad = find(isnan(sortedPDI));
sortedNP(bad) = [];
sortedPDI(bad) = [];


figure,
plot(sortedNP, (1:length(sortedNP))', 'ro')
hold on
plot(sortedPDI, 1:length(sortedPDI), 'k.')
title(region)
ylabel('Trial Number')
xlabel('Time')
%legend({'Nosepoke Exit', 'Significant PV divergence'})

figtitle = ['PVDivergence-NPlength_sessions_raster_',region];
figfile = [figDir,'PopSelectivity\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);


% figure,
% plot([zeros(length(allNPoff),1), sortedNP], [(1:length(allNPoff))', (1:length(allNPoff))'], 'ro')
% hold on
% plot(sortedPeaks, 1:length(allPeakTimes), 'k.')
% title(region)
% ylabel('Trial Number')
% 
% figtitle = 'PeakPDI-NPlength_raster';
% figfile = [figDir,'PopSelectivity\',figtitle];
% %saveas(gcf,figfile,'fig');
% print('-dpdf', figfile);
% print('-djpeg', figfile);
% 
% % figure
% % plot([zeros(length(allNPoff),1), allNPoff], [(1:length(allNPoff))', (1:length(allNPoff))'], 'rx')
% % hold on
% % plot(allDivTimes, 1:length(allDivTimes), 'go')

%% FIgure 2 - correlation

% keepinds = find(~isnan(allDivTimes));
% 
% allDivTimes = allDivTimes(keepinds);
% allNPoff = allNPoff(keepinds);

figure
plot(sortedNP, sortedPDI, 'k.')
hold on

%keepinds = ~isnan(allDivTimes);
fit = polyfit(sortedPDI,sortedNP,1);
plot([min(sortedPDI), max(sortedPDI)], polyval(fit,[min(sortedPDI), max(sortedPDI)] ))

xlabel('First significant divergence time')
ylabel('Nosepoke exit time')
title(region)
[CCdiv,p_div] = corrcoef(sortedPDI,sortedNP);
R = CCdiv(1,2);
p = p_div(1,2);

text(0.8,.6,['R = ',num2str(R), newline, 'p = ', num2str(p)])

figtitle = ['PVDivergence-NPlength_sessions_correlation_',region];
figfile = [figDir,'PopSelectivity\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);


% figure
% plot(allPeakTimes, allNPoff, 'k.')
% fit = polyfit(allPeakTimes,allNPoff,1);
% hold on
% plot([min(allPeakTimes), max(allPeakTimes)], polyval(fit,[min(allPeakTimes), max(allPeakTimes)] ))
% xlabel('Peak divergence time')
% ylabel('Time before nosepoke exit')
% title(region)
% 
% [CCpeak, p_peak] = corrcoef(allPeakTimes,allNPoff)
% 
% figtitle = 'PeakPDI-NPlength_correlation';
% figfile = [figDir,'PopSelectivity\',figtitle];
% %saveas(gcf,figfile,'fig');
% print('-dpdf', figfile);
% print('-djpeg', figfile);
end