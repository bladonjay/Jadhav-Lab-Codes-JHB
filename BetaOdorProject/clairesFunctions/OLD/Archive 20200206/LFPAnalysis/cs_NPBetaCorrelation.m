%raster showing timing of odor onset, neural activity divergence (PDI), and
%nosepoke exit

%fig1: raster
%fig2: correlation between PDI and NP exit
%fig3: comparison between CA1 and PFC mean digergence time, and time
%between divergence and NP exit (with STD) - may show which one is more
%reliable?
clear
close all
[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};
%animals = {'CS33'};
regions = {'CA1','PFC'};

for r = 1:length(regions)
    region = regions{r};
    
    allBeta = [];
    allNPoff = [];
    allDivTimes = [];
    %allPeakTimes = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animaldir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        files = dir([animaldir,animal,'trialSigPDI_',region,'*']);
        
        load([animaldir,animal,'betaWindows.mat']);
        
        if ~isempty(files)
            for f = 1:length(files)
                load(files(f).name);
                day = length(trialSigPDI);
                daystr = getTwoDigitNumber(day);
                load([animaldir,animal,'nosepokeWindow',daystr,'.mat'])
                load([animaldir,animal,'odorTriggers',daystr,'.mat'])
                
                %load([animaldir,animal,'trialSigPDI_',region,daystr,'.mat'])

                epochs = find(~cellfun(@isempty,trialSigPDI{day}));
                for e = 1:length(epochs)
                    epoch = epochs(e);
                    
                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                    correctinds = sort([correct_left;correct_right]);
                    
                    windows = nosepokeWindow{day}{epoch}(correctinds,:);
                    
                    npoff = windows(:,2)-windows(:,1); %normalize time windows
                    divTimes = trialSigPDI{day}{epoch}(correctinds,:);
                    %peakTimes = trialPeakPDI{day}{epoch}(correctinds,:);
                    betaOnsets = betaWindows{day}{epoch}(:,1);
                    
                    betaOnsets = betaOnsets(correctinds) - windows(:,1); %normalize
                    
                    allBeta = [allBeta; betaOnsets];
                    allNPoff = [allNPoff; npoff];
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

betainds = find(allBeta < 1.5);
allBeta = allBeta(betainds);
allDivTimes = allDivTimes(betainds);

pvinds = find(allDivTimes < 1.5);
allDivTimes = allDivTimes(pvinds);
allBeta = allBeta(pvinds);


[sortedPDI, ind] = sort(allDivTimes);
% sortedPDI = allDivTimes(ind);
sortedBeta = allBeta(ind);
%sortedPeaks = allPeakTimes(ind);

keepinds = ~isnan(sortedPDI);
sortedPDI = sortedPDI(keepinds);
%sortedNP = sortedNP(keepinds);
sortedBeta = sortedBeta(keepinds);


bad = find(sortedBeta < 0.01); %exclude trials where beta started off high
sortedBeta(bad) = nan;
keepinds2 = ~isnan(sortedBeta);
sortedBeta = sortedBeta(keepinds2);
sortedPDI = sortedPDI(keepinds2);
%sortedNP = sortedNP(keepinds2);





% figure,
% %plot([zeros(length(sortedNP),1), sortedNP], [(1:length(sortedNP))', (1:length(sortedNP))'], 'ro')
% plot( sortedNP, (1:length(sortedNP))', 'ro')
% hold on
% plot(sortedPDI, 1:length(sortedPDI), 'k.')
% plot(sortedBeta, 1:length(sortedBeta), 'gx')
% title(region)
% ylabel('Trial Number')
% xlabel('Time')
% %legend({'Nosepoke Exit', 'Significant PV divergence'})
% 
% figtitle = ['PVDivergence-NPlength_raster_',region];
% figfile = [figDir,'PopSelectivity\',figtitle];
% %saveas(gcf,figfile,'fig');
% print('-dpdf', figfile);
% print('-djpeg', figfile);


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

%% FIgure 2 - PV divergence time/ Beta Onset Time correlation

% keepinds = find(~isnan(allDivTimes));
% 
% allDivTimes = allDivTimes(keepinds);
% allNPoff = allNPoff(keepinds);

figure
plot(sortedPDI, sortedBeta, 'k.')
hold on

%fit = polyfit(allDivTimes(keepinds),allNPoff(keepinds),1);
fit = polyfit(sortedPDI,sortedBeta,1);
plot([min(sortedPDI), max(sortedPDI)], polyval(fit,[min(sortedPDI), max(sortedPDI)] ))

xlabel('First significant divergence time')
ylabel('Beta Onset Time')
title(region)
[CCdiv,p_div] = corrcoef(sortedPDI,sortedBeta);
R = CCdiv(1,2);
p = p_div(1,2);

text(1,1,['R = ',num2str(R), newline, 'p = ', num2str(p)])

figtitle = ['PVDivergence-BetaOnset_correlation_',region];
figfile = [figDir,'PopSelectivity\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);



%% Figre 3 - Beta Onset/ NP length correlation

% figure
% plot(sortedBeta, sortedNP, 'k.')
% hold on
% 
% %fit = polyfit(allDivTimes(keepinds),allNPoff(keepinds),1);
% fit = polyfit(sortedBeta,sortedNP, 1);
% plot([min(sortedBeta), max(sortedBeta)], polyval(fit,[min(sortedBeta), max(sortedBeta)] ))
% 
% xlabel('Beta Onset Time')
% ylabel('Nosepoke Disnegagement Time')
% title(region)
% [CCdiv,p_div] = corrcoef(sortedBeta, sortedNP);
% R = CCdiv(1,2);
% p = p_div(1,2);
% 
% text(1,1.4,['R = ',num2str(R), newline, 'p = ', num2str(p)])
% 
% figtitle = ['BetaOnset-NPtime_correlation_',region];
% figfile = [figDir,'PopSelectivity\',figtitle];
% %saveas(gcf,figfile,'fig');
% print('-dpdf', figfile);
% print('-djpeg', figfile);

end