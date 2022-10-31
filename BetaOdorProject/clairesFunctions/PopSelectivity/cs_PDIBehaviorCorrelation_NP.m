%raster showing timing of odor onset, neural activity divergence (PDI), and
%nosepoke exit

%fig1: raster
%fig2: correlation between PDI and NP exit
%fig3: comparison between CA1 and PFC mean digergence time, and time
%between divergence and NP exit (with STD) - may show which one is more
%reliable?
close all
clear
[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
regions = {'CA1','PFC'};

for r = 1:length(regions)
    region = regions{r};
    
    allNPoff = [];
    allDivTimes = [];
    allPeakTimes = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animaldir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        files = dir([animaldir,animal,'trialSigPDI_NP_',region,'*']);
        
        if ~isempty(files)
            for f = 1:length(files)
                load([animaldir,files(f).name]);
                day = length(trialSigPDI);
                daystr = getTwoDigitNumber(day);
                load([animaldir,animal,'nosepokeWindow',daystr,'.mat'])
                load([animaldir,animal,'odorTriggers',daystr,'.mat'])
                
                %load([animaldir,animal,'trialSigPDI_',region,daystr,'.mat'])

                epochs = find(~cellfun(@isempty,trialSigPDI{day}));
                for e = 1:length(epochs)
                    epoch = epochs(e);
                    
                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                    correctinds = [correct_left;correct_right];
                    
                    
                    windows = nosepokeWindow{day}{epoch}(correctinds,:);
                    
                    npoff = windows(:,2)-windows(:,1); %normalize time windows
                    divTimes = trialSigPDI{day}{epoch};
                    
                    %remove any trials where np wasn't long enough
                    test = npoff<0.5;
                    if sum(test) >0
                       npoff = npoff(~test);
                       divTimes = divTimes(~test);
                    end
                    
                    %peakTimes = trialPeakPDI{day}{epoch}(correctinds,:);
                    
                    allNPoff = [allNPoff; npoff];
                    allDivTimes = [allDivTimes; divTimes];
                    %allPeakTimes = [allPeakTimes; peakTimes];
                end
            end
        end
    end
%% Figure 1 - raster
 npinds = find(allNPoff < 1);
 allNPoff = allNPoff(npinds);
 allDivTimes = allDivTimes(npinds);

% pvinds = find(allDivTimes < 1.5);
% allDivTimes = allDivTimes(pvinds);
% allNPoff = allNPoff(pvinds);


[sortedNP, ind] = sort(allNPoff);
sortedPDI = allDivTimes(ind);
%sortedPeaks = allPeakTimes(ind);

keepinds = ~isnan(sortedPDI);
sortedPDI = sortedPDI(keepinds);
sortedNP = sortedNP(keepinds);

figure, 
%plot([zeros(length(sortedNP),1), sortedNP], [(1:length(sortedNP))', (1:length(sortedNP))'], 'ro')
barh(1:length(sortedNP),sortedNP,1,'FaceColor','r','FaceAlpha',0.5);
%plot( sortedNP, (1:length(sortedNP))', 'r.','MarkerSize',12)
hold on
box off
plot(sortedPDI, 1:length(sortedPDI), 'k.','MarkerSize',12)
title(region)
axis([0, 1, 0, length(sortedNP)+5])
%set(gcf,'Position',[460, 650, 420, 420])
ylabel('Trial Number')
xlabel('Time')
%legend({'Nosepoke Exit', 'Significant PV divergence'})

figtitle = ['PVDivergence-NPlength_raster_',region];
figfile = [figDir,'PopSelectivity\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);



%% FIgure 2 - correlation
% 
% figure, 
% plot(sortedPDI, sortedNP, 'k.','MarkerSize',12)
% box off
% hold on
% 
% %fit = polyfit(allDivTimes(keepinds),allNPoff(keepinds),1);
% fit = polyfit(sortedPDI,sortedNP,1);
% plot([min(sortedPDI), max(sortedPDI)], polyval(fit,[min(sortedPDI), max(sortedPDI)] ))
% 
% set(gcf,'Position',[1200,300,510,420]);
% xlabel('First significant divergence time')
% ylabel('Nosepoke exit time')
% title(region)
[CCdiv,p_div] = corrcoef(sortedPDI,sortedNP);
R = CCdiv(1,2);
p = p_div(1,2);

% text(0.2,1,['R = ',num2str(R), newline, 'p = ', num2str(p)])
% 
% figtitle = ['PVDivergence-NPlength_correlation_',region];
% figfile = [figDir,'PopSelectivity\',figtitle];
% %saveas(gcf,figfile,'fig');
% print('-dpdf', figfile);
% print('-djpeg', figfile);



%sextiles
figure
[means] = cs_sextiles(sortedNP,sortedPDI);

 ylabel('PV divergence')
 xlabel('Decision Latency sextile')
 text(0,means(1),['R = ',num2str(R), 'p = ',num2str(round(p,2,'significant'))]);
 figtitle = ['PVDivergence-NPlength_sextile_',region];
figfile = [figDir,'PopSelectivity\',figtitle];
xlim([0 7])
%ylim([0.4 0.7])
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);

end