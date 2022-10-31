%raster showing timing of odor onset, neural activity divergence (PDI), and
%nosepoke exit

%fig1: raster
%fig2: correlation between PDI and NP exit
%fig3: comparison between CA1 and PFC mean digergence time, and time
%between divergence and NP exit (with STD) - may show which one is more
%reliable?
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};
%animals = {'CS33'};
regions = {'CA1','PFC'};
figure, hold on
%for r = 1:length(regions)
    %region = regions{r};
    
 
    allDivTimes = [];

    for a = 1:length(animals)
        animal = animals{a};
        
        animaldir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        PFCfiles = dir([animaldir,animal,'trialSigPDI_PFC*']);
        CA1files = dir([animaldir,animal,'trialSigPDI_CA1*']);
        
        PFCfiles = {PFCfiles.name};
        CA1files = {CA1files.name};
        %combine = {PFCfiles.name,CA1files.name};
        PFCdays = cellfun(@(x) str2double(x(21)),PFCfiles);
        CA1days = cellfun(@(x) str2double(x(21)),CA1files);
        
        days = intersect(PFCdays,CA1days);
        
        
        
        if ~isempty(days)
            for f = 1:length(days)
                day = length(days(d));
                daystr = getTwoDigitNumber(day);
                
                load(['trialSigPDI_PFC',daystr,'.mat']);

                load([animaldir,animal,'odorTriggers',daystr,'.mat'])

                epochs = find(~cellfun(@isempty,trialSigPDI{day}));
                for e = 1:length(epochs)
                    epoch = epochs(e);
                    
                    [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                    correctinds = sort([correct_left;correct_right]);
                    
                    divTimes = trialSigPDI{day}{epoch}(correctinds,:);
                
                    
                    
                    allDivTimes = [allDivTimes; divTimes];

                end
            end
        end
    end
%% Figure 1 - Histogram

%allDivTimes = allDivTimes(~isnan(allDivTimes));


alldiv.(region) = allDivTimes;
%histogram(allDivTimes,30)

%mn = mean(allDivTimes);
%plot([mn, mn], [0 30], '--');

%end
CA1 = alldiv.CA1;
PFC = alldiv.PFC;

[p] = ranksum(CA1, PFC,0.05)
[~,p] = ttest2(CA1, PFC,0.05)

[CA1,ind] = sort(CA1);
PFC = interp1(1:numel(PFC),PFC,linspace(1,numel(PFC),numel(CA1)));
[R,p] = corrcoef(CA1,PFC)


% 
% figtitle = 'PeakPDI-NPlength_correlation';
% figfile = [figDir,'PopSelectivity\',figtitle];
% %saveas(gcf,figfile,'fig');
% print('-dpdf', figfile);
% print('-djpeg', figfile);