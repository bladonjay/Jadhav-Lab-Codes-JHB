% Replications PV differential index calculation from Igarashi et al 2012.

% win = [0.2 1.2];
% binsize = 0.1;


function cs_individualTrialPDI_peak(win,binsize)
close all
%% Params
animals = {'CS31','CS33','CS34','CS35'};

[topDir] = cs_setPaths();
regions = {'CA1','PFC'};



winstring = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
binstr = [num2str(binsize),'msBins'];
bins = (-win(1):binsize:win(2));

load([topDir, 'AnalysesAcrossAnimals\significantPDI.mat'])

%% 
for r = 1:length(regions)
    region = regions{r};
    
    %load selective cells only
    load([topDir, 'AnalysesAcrossAnimals\selectiveCells_', region, '.mat'])
    uniquedays = unique(selectivecells(:,1:2),'rows');
    
    %use value calculated from cs_PDI (from all cells/animals)
    load([topDir, 'AnalysesAcrossAnimals\significantPDI.mat'])
    sigPDI = significantPDI.(region)(2);
    
    for ad = 1:size(uniquedays,1)
        
        cells = selectivecells(ismember(selectivecells(:,1:2), uniquedays(ad,:), 'rows'),:);
        if size(cells,1) > 2
            
            animal = animals{uniquedays(ad,1)};
            
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            
            day = uniquedays(ad,2);
            daystr = getTwoDigitNumber(day);
            
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            load([animDir,animal,'spikes',daystr,'.mat'])
            
            runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
            for ep = 1:length(runeps)
                
                sessionSpikingData = {};
                epoch = runeps(ep);
                epstr = getTwoDigitNumber(epoch);
                disp(['Doing animal ', animal, ' day ', daystr, ' epoch ', epstr]);
                
                [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                
                trigs = odorTriggers{day}{epoch}.allTriggers;
                
                lefttrigs = trigs(correct_left);
                righttrigs = trigs(correct_right);
                
                
               
%% Calculate the session PDI using the means

                
                for tr = 1:length(trigs)
                    trig = trigs(tr);
                    
                    trigwin = [trig-win(1), trig + win(2)];
                    trigbins = [trig-win(1):binsize:trig + win(2)];
      
                    %create one large cell array, then use later when
                    %looking at individual trials
                    for c = 1:size(cells,1)
                        if ismember(trig, lefttrigs) || ismember(trig, righttrigs) %only correct trials
                            
                            if ~isempty(spikes{day}{epoch}{cells(c,3)}{cells(c,4)}) 
                                %gather spikes from each trial
                                cellspikes = spikes{day}{epoch}{cells(c,3)}{cells(c,4)}.data(:,1);
                                winspikes = cellspikes(cellspikes > trigwin(1) & cellspikes <= trigwin(2));
                                binspikecount = histcounts(winspikes',trigbins);
                                trialFR = binspikecount./binsize;
                                trialFR = smoothdata(trialFR, 'gaussian', 5);
                                
                                sessionSpikingData{c}(tr,:) = trialFR;
                            else
                                sessionSpikingData{c}(tr,:) = zeros(1, length(bins)-1);
                            end
                            
                        else
                            sessionSpikingData{c}(tr,:) = nan(1, length(bins)-1);
                        end
                        
                        
                    end
                    
                end
                
                %separate into left and right trials
                leftSpikingData = cellfun(@(x) x(correct_left,:), sessionSpikingData, 'UniformOutput', false);
                rightSpikingData = cellfun(@(x) x(correct_right,:), sessionSpikingData, 'UniformOutput', false);
                
                %get mean PVs across trials
                leftmeanPVs = vertcat(cell2mat(cellfun(@(x) nanmean(x, 1), leftSpikingData,'UniformOutput',false)'));
                rightmeanPVs = vertcat(cell2mat(cellfun(@(x) nanmean(x, 1), rightSpikingData,'UniformOutput',false)'));
                
                
%% compare individual trials to the session means
                for tr = 1:length(trigs)
                    trig = trigs(tr);
                    %timepoint = 0;
                    if ismember(trig, lefttrigs) || ismember(trig, righttrigs) %only correct trials
                        trialPVs = vertcat(cell2mat(cellfun(@(x) x(tr,:), sessionSpikingData, 'UniformOutput',false)'));
                        
                        %compare to the mean of the opposite side trials
                        if ismember(trig, lefttrigs)
                            oppositePVs = rightmeanPVs;
                        elseif ismember(trig, righttrigs)
                            oppositePVs = leftmeanPVs;
                        end
                        
                        %calc trial PDI
                        trialPDI = zeros(1,length(bins)-1);
                        for b = 1:length(bins)-1
                            trialPV = trialPVs(:,b);
                            oppositePV = oppositePVs(:,b);
                            
                            CCbin = corrcoef(trialPV, oppositePV);
%                              mag1 = sqrt(sum(trialPV.^2));
%                             mag2 = sqrt(sum(oppositePV.^2));
%                             
%                             CCbin = (dot(trialPV,oppositePV))/(mag1 * mag2);
                            
                            PDIbin = 1-abs(CCbin(1,2)); %PV Differential Index, from Igarashi et al 2014
                            if isnan(PDIbin)
                                PDIbin = 0;
                            end
                            
                            trialPDI(b) = PDIbin;
                            
                        end
                        %plot(trialPDI)
                        %pause,
%% determine time of PEAK PDI
                        binsaftertrig = find(bins(1:end-1)> 0.1); %at least 200ms after NP
                        [maxpdi, maxbin] = max(trialPDI(binsaftertrig));

                            if maxbin == length(binsaftertrig)%if its the last bin, use the second highest
                                tmp = trialPDI(binsaftertrig);
                                tmp(maxbin) = nan;
                                [maxpdi, maxbin] = max(tmp);
                            end
                            
                            
                            if maxpdi >= sigPDI %max must be higher than sigPDI
                            
                            maxbin = maxbin  + (binsaftertrig(1)-1);
                            prevbin = maxbin-1;
                            nextbin = maxbin+1;
                            
                            binstouse = [prevbin, maxbin, nextbin];
                            
                            fit = polyfit(bins(binstouse), trialPDI(binstouse), 2);
                            
                            times = linspace(bins(prevbin), bins(nextbin));
                            
                            
                            
                            [peakPDI, ind] = max(polyval(fit,times));
                            
                            timepoint = times(ind);
%                             
                            
                            
                            
%                             disp(timepoint)
%                             figure, plot(bins(1:end-1), trialPDI)
%                             hold on
%                             plot([bins(1) bins(end-1)], [sigPDI sigPDI], 'k--');
%                             
%                             fitpnts = polyval(fit,times);
%                             plot(times, fitpnts, 'b--');
%                             
%                             plot(timepoint, sigPDI, 'rx');
%                             pause; close all
% %                             
                            trialPeakPDI{day}{epoch}(tr,1) = timepoint;
                        else
                            trialPeakPDI{day}{epoch}(tr,1) = NaN;
                        end
                        
                    else %incorrect trials
                        trialPeakPDI{day}{epoch}(tr,1) = NaN;
                    end
                end
            end
%% Save
            save([animDir,animal,'trialPeakPDI_',region,daystr,'.mat'],'trialPeakPDI')
            clear trialPeakPDI
        end
        
    end
end
end