% Replications PV differential index calculation from Igarashi et al 2012.
clear
win = [0 2];
binsize = 0.1;

trialtype = 'correct';
% mintrialspikes = 50; %minimum spikes across all cells on trial


g = 7; %smoothing %use 7
winsize = 2; %2
%function cs_individualTrialPDI(win,binsize)
close all
%% Params
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir,figDir] = cs_setPaths();


regions = {'CA1','PFC'};



winstring = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
binstr = [num2str(binsize),'msBins'];
bins = (-win(1):binsize:win(2));

load([topDir, 'AnalysesAcrossAnimals\significantPDI.mat'])

%%
for r = 1:length(regions)
    region = regions{r};
    switch trialtype
        case 'correct'
            trialstr = '';
            switch region
                case 'CA1'
                    mintrialspikes = 40;
                case 'PFC'
                    mintrialspikes = 40;
            end
        case 'incorrect'
            trialstr = 'incorrect_';
            switch region
                case 'CA1'
                    mintrialspikes = 50;
                case 'PFC'
                    mintrialspikes = 20;
            end
    end
    vals = [];
    %delete previous files
%     for a = 1:length(animals)
%         animal = animals{a};
%         animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
%         files = dir([animDir,animal,'trialSigPDI_',trialstr,region,'*']);
%         for f = 1:length(files)
%             filename = files(f).name;
%             delete([animDir,filename])
%         end
%     end
    
    %load selective cells only
    load([topDir, 'AnalysesAcrossAnimals\selectiveCells_', region, '.mat'])
    uniquedays = unique(selectivecells(:,1:2),'rows');
    
    %use value calculated from cs_PDI (from all cells/animals)
    load([topDir, 'AnalysesAcrossAnimals\significantPDI.mat'])
    sigPDI = significantPDI.(region)(2);
    
    for ad = 1:size(uniquedays,1)
        
        cells = selectivecells(ismember(selectivecells(:,1:2), uniquedays(ad,:), 'rows'),:);
        if size(cells,1) > 2 %need at least 3 cells to calc PDI
            
            animal = animals{uniquedays(ad,1)};
            
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            
            day = uniquedays(ad,2);
            daystr = getTwoDigitNumber(day);
            
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            load([animDir,animal,'nosepokeWindow',daystr,'.mat'])
            load([animDir,animal,'spikes',daystr,'.mat'])
            
            runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
            for ep = 1:length(runeps)
                
                
                sessionSpikingData = {};
                epoch = runeps(ep);
                epstr = getTwoDigitNumber(epoch);
                disp(['Doing animal ', animal, ' day ', daystr, ' epoch ', epstr]);
                
                [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                
                npwin = nosepokeWindow{day}{epoch};
                
                switch trialtype
                    case 'correct'
                        lefttrigs = npwin(correct_left,:);
                        righttrigs = npwin(correct_right,:);
                    case 'incorrect'
                        lefttrigs = npwin(incorrect_left,:);
                        righttrigs = npwin(incorrect_right,:);
                end
                
                
                %% Calculate the session PDI using the means
                
                
                for n = 1:length(npwin)
                    npon = npwin(n,1);
                    
                    trigwin = [npon-win(1), npon + win(2)];
                    trigbins = [npon-win(1):binsize:npon + win(2)];
                    
                    %create one large cell array, then use later when
                    %looking at individual trials
                    for c = 1:size(cells,1)
                        if ismember(npon, lefttrigs) || ismember(npon, righttrigs) %only correct trials
                            
                            if ~isempty(spikes{day}{epoch}{cells(c,3)}{cells(c,4)}.data)
                                %gather spikes from each trial
                                cellspikes = spikes{day}{epoch}{cells(c,3)}{cells(c,4)}.data(:,1);
                                winspikes = cellspikes(cellspikes > trigwin(1) & cellspikes <= trigwin(2));
                                %spiking criterion for each trial- must be
                                %more than 5 spikes
                                binspikecount = histcounts(winspikes',trigbins);
                                
                                %if sum(binspikecount) > 3
                                %                                     trialFR = binspikecount./binsize;
                                %s = gaussian(1,5);
                                %trialFR = filter(s,1,trialFR);
                                %plot(trialFR);
                                %                                     trialFR = smoothdata(trialFR, 'gaussian', 8);
                                %sessionSpikingData{c}(tr,:) = trialFR;
                                %                                     plot(trialFR);
                                sessionSpikingData{c}(n,:) = binspikecount;
                                %else
                                %sessionSpikingData{c}(tr,:) = nan(1, length(bins)-1);
                                %end
                                
                            else
                                sessionSpikingData{c}(n,:) = nan(1, length(bins)-1);
                            end
                            
                        else
                            sessionSpikingData{c}(n,:) = nan(1, length(bins)-1);
                        end
                        
                        
                    end
                    
                end
                
                for d = 1:size(sessionSpikingData{1},1)
                    trialspikes = nansum(cell2mat(cellfun(@(x) nansum(x(d,:)), sessionSpikingData, 'UniformOutput', false)));
                    
                    if trialspikes < mintrialspikes %minimum spikes from all cells on trial
                        
                        for c = 1:length(sessionSpikingData)
                            
                            sessionSpikingData{c}(d,:) = nan;
                        end
                    end
                    
                end
                
                %trialFR = binspikecount./binsize;
                
                %separate into left and right trials
                switch trialtype
                    case 'correct'
                        leftSpikingData = cellfun(@(x) x(correct_left,:), sessionSpikingData, 'UniformOutput', false);
                        rightSpikingData = cellfun(@(x) x(correct_right,:), sessionSpikingData, 'UniformOutput', false);
                    case 'incorrect'
                        leftSpikingData = cellfun(@(x) x(incorrect_left,:), sessionSpikingData, 'UniformOutput', false);
                        rightSpikingData = cellfun(@(x) x(incorrect_right,:), sessionSpikingData, 'UniformOutput', false);
                end
                %get mean PVs across trials
                
                
                leftmeanPVs = vertcat(cell2mat(cellfun(@(x) nanmean(x, 1), leftSpikingData,'UniformOutput',false)'));
                rightmeanPVs = vertcat(cell2mat(cellfun(@(x) nanmean(x, 1), rightSpikingData,'UniformOutput',false)'));
                
                leftmeanPVs = smoothdata(leftmeanPVs./binsize,2,'gaussian',g);
                rightmeanPVs = smoothdata(rightmeanPVs./binsize,2,'gaussian',g);
                % %
                %
                %if any rows are all NaN, then those cells never crossed
                %spiking threshold. Remove these rows so they are not
                %considered in PV analysis
                %                                 leftmeanPVs = leftmeanPVs(any(leftmeanPVs,2),:);
                %                                 rightmeanPVs = rightmeanPVs(any(rightmeanPVs,2),:);
                %
                %                                 if size(leftmeanPVs,1) < 2 | size(rightmeanPVs,1) <2
                %                                     flag = 1;
                %                                     continue
                %                                 end
                
                %% compare individual trials to the session means
                for n = 1:length(npwin)
                    npon = npwin(n,1);
                    npoff = npwin(n,2)-npon;
                    %timepoint = 0;
                    if ismember(npon, lefttrigs) || ismember(npon, righttrigs) %only correct trials
                        trialspiking = vertcat(cell2mat(cellfun(@(x) x(n,:), sessionSpikingData, 'UniformOutput',false)'));
                        trialPVs = trialspiking./binsize;
                        trialPVs = smoothdata(trialPVs,2, 'gaussian', g);
                        
                        if all(any(isnan(trialPVs),2)) || all(any(isnan(leftmeanPVs),2)) || all(any(isnan(rightmeanPVs),2))
                            trialSigPDI{day}{epoch}(n,1) = NaN;
                            continue
                        else
                            goodtrial = find(~any(isnan(trialPVs),2));
                            goodmnleft = find(~any(isnan(rightmeanPVs),2));
                            goodmnright = find(~any(isnan(leftmeanPVs),2));
                            
                            mninds = intersect(goodmnleft, goodmnright);
                            rowinds = intersect(mninds, goodtrial);
                            %rowinds = find(~any(isnan(trialPVs),2));
                            
                            trialPVs = trialPVs(rowinds,:);
                            if length(rowinds) < 2
                               continue
                            end
                        end
                        %testSpikingCriterion = any(trialPVs,2);
                        %if sum(testSpikingCriterion) < 2
                        %    trialSigPDI{day}{epoch}(tr,1) = NaN;
                        %    continue
                        %else
                        %    trialPVs = trialPVs(testSpikingCriterion,:);
                        %end
                        
                        %compare to the mean of the opposite side trials
                        if ismember(npon, lefttrigs)
                            %oppositePVs = rightmeanPVs(testSpikingCriterion,:);
                            oppositePVs = rightmeanPVs(rowinds,:);
                        elseif ismember(npon, righttrigs)
                            %oppositePVs = leftmeanPVs(testSpikingCriterion,:);
                            oppositePVs = leftmeanPVs(rowinds,:);
                        end
                        
                        %calc trial PDI
                        trialPDI = zeros(1,length(bins)-1);
                        for b = 1:length(bins)-1
                            trialPV = trialPVs(:,b);
                            oppositePV = oppositePVs(:,b);
                            
                            CCbin = corrcoef(trialPV, oppositePV);
                            %
                            
                            PDIbin = 1-abs(CCbin(1,2)); %PV Differential Index, from Igarashi et al 2014
                            if isnan(PDIbin)
                                PDIbin = 0;
                            end
                            
                            trialPDI(b) = PDIbin;
                            
                        end
                        
                        %% Get Bin During Which He Exited NP
                        [ind] = lookup(npoff,bins);
                        npexitbin = bins(ind);
                        if npoff < npexitbin
                            ind1 = ind-1;
                            ind2 = ind;
                        elseif npoff > npexitbin
                            ind1 = ind;
                            ind2 = ind+1;
                        end
                        
                        bin1 = bins(ind1);
                        bin2 = bins(ind2);
                        
                        val = interp1([bin1, bin2], [trialPDI(ind1), trialPDI(ind2)], npoff);
                        
                        if ~isnan(val) & val > 0.01
                            vals = [vals; val];
                            if val < 0.1
%                                 plot(bins(1:end-1), trialPDI);
%                                 hold on
%                                 plot([npoff npoff],[0 1],'r--')
%                                 text(0.2, 0.7, [region,' ',num2str(val)]);
%                                 keyboard
%                                 close all
                            end
                        end
                        
                        
                        
                        
                        %                     else %incorrect trials
                        %                         trialSigPDI{day}{epoch}(n,1) = NaN;
                    end
                end
            end
            %% Save
%             save([animDir,animal,'trialSigPDI_',trialstr,region,daystr,'.mat'],'trialSigPDI')
%             clear trialSigPDI
        end
        
    end
    allvals.(region) = vals;
end

mnCA1 = mean(allvals.CA1);
mnPFC = mean(allvals.PFC);

CA1 = allvals.CA1;
PFC = allvals.PFC;

figure, hold on
bplot(CA1,1,'nooutliers','nomean','width',0.5);
bplot(PFC,2,'nooutliers','nomean','width',0.5);
xticks([1 2])
xticklabels(regions);
axis([0 3 0 1])

[p] = ranksum(allvals.CA1, allvals.PFC,0.05);

text(2,0.9,['p = ',num2str(p)]);
figtitle = 'PVDivergenceValueComparison';
figfile = [figDir,'PopSelectivity\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);


% 
% figure, hold on
% scatter(ones(length(allvals.CA1),1), allvals.CA1, 15, 'filled', 'jitter', 'on', 'jitterAmount', 0.05);
% scatter([ones(length(allvals.PFC),1)+1], allvals.PFC, 15, 'filled', 'jitter', 'on', 'jitterAmount', 0.05);
% axis([0.5 2.5 0 1.1])
% 
% figure, bar([1 2], [mnCA1,mnPFC]);
% semCA1 = std(CA1)/sqrt(length(CA1));   
% semPFC = std(PFC)/sqrt(length(PFC));
% hold on
% er = errorbar([1,2],[mnCA1, mnPFC],[mnCA1-semCA1, mnPFC-semPFC],[mnCA1+semCA1, mnPFC+semPFC]);


% figure, hold on
% [X,Y] = cumhist(CA1, [min(CA1) max(CA1)], 1);
% [X,Y] = cumhist(PFC, [min(PFC) max(PFC)], 1);
% 
% figure, hold on
% histogram(allvals.CA1,30)
% histogram(allvals.PFC,30)