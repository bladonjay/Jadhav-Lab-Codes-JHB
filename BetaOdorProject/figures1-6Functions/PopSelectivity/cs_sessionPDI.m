% Replications PV differential index calculation from Igarashi et al 2012.

win = [0 1.5];
binsize = 0.1;


%cs_session(win,binsize)
close all
%% Params
animals = {'CS31','CS33','CS34','CS35'};

[topDir] = cs_setPaths();
regions = {'CA1'};



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
                
                
               
%% Calculate the session PDI
          
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
                                 s = gaussian(3,6);
                                      trialFR = filter(s,1,trialFR); 
                                %trialFR = smoothdata(trialFR, 'gaussian', 5);
                                
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
                
                
                for b = 1:length(bins)-1
                    PV_l = leftmeanPVs(:,b);
                    PV_r = rightmeanPVs(:,b);
                    CCbin = corrcoef(PV_l, PV_r);
                    PDIbin = 1-abs(CCbin(1,2));
                    
                    if isnan(PDIbin)
                        PDIbin = 0;
                    end
                    
                    PDI(b) = PDIbin;
                end
                
                figure, plot(PDI);
                
                %% Determine where session PDI reaches significance threshold
                binsaftertrig = find(bins(1:end-1)> 0); %at least 200ms after NP
                sigbin = find(PDI(binsaftertrig) >= sigPDI, 1, 'first');
                
                if ~isempty(sigbin) %if it becomes significant at some point during the trial
                             
                            sigbin = sigbin  + (binsaftertrig(1)-1);
                            
                            sigval = PDI(sigbin);
                            prevval = PDI(sigbin-1);
                            t1 = bins(sigbin-1);
                            t2 = bins(sigbin);
                            
                            %interpolate to find the exact crossing
                            %time where the PDI becomes significant
                            timepoint = interp1([prevval, sigval], [t1, t2], sigPDI); 
                            
                            while (sigval < prevval || sigval <  sigPDI || prevval > sigPDI)
                                    sigbin = sigbin+1;
                                if sigbin <= length(PDI)
                                    sigval = PDI(sigbin);
                                    prevval = PDI(sigbin-1);
                                    t1 = bins(sigbin-1);
                                    t2 = bins(sigbin);
                                    if ~isequal(prevval, sigval)
                                        timepoint = interp1([prevval, sigval], [t1, t2], sigPDI); 
                                    end
                                else 
                                    timepoint = nan;
                                    break %stops the loop
                                end
                            end
                                sessionSigPDI{day}{epoch} = timepoint;
                        else
                            sessionSigPDI{day}{epoch} = NaN;
                        end
            end
                %% Save
                save([animDir,animal,'sessionSigPDI_',region,daystr,'.mat'],'sessionSigPDI')
                clear sessionSigPDI
        end
            
        end
        
end