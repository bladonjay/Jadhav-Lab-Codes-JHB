% Replications PV differential index calculation from Igarashi et al 2012.
%changed the way smoothing is done- use interpolation 
clear
close all
%% Params
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir] = cs_setPaths();
regions = {'CA1','PFC'};
win = [0 1.0];
binsize = 0.1;
interpfactor = 5;
celltype = 'np'; %np or selective
trialtype = 'correct';
minensemblecells = 5;
% mintrialspikes = 50; %minimum spikes across all cells on trial

bins = (-win(1):binsize:win(2));
newbins = bins(1):binsize/interpfactor:bins(end)-binsize/interpfactor;

winstring = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
binstr = [num2str(binsize),'msBins'];

load([topDir, 'AnalysesAcrossAnimals\significantPDI_NP.mat'])

%%
for r = 1:length(regions)
    region = regions{r};
    switch trialtype
        case 'correct'
            trialstr = '';
            switch region
                case 'CA1'
                    mintrialspikes = 10;
                case 'PFC'
                    mintrialspikes = 10;
            end
        case 'incorrect'
            trialstr = 'incorrect_';
            switch region
                case 'CA1'
                    mintrialspikes = 10;
                case 'PFC'
                    mintrialspikes = 10;
            end
    end
   vals = [];
    %delete previous files
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    files = dir([animDir,animal,'trialSigPDI_NP_',trialstr,region,'*']);
    for f = 1:length(files)
        filename = files(f).name;
        delete([animDir,filename])
    end
end
    
    %load cells
    load([topDir, 'AnalysesAcrossAnimals\',celltype,'Cells_', region, '.mat'])
    
    try
        eval(['allcells = ',celltype,'Cells;']);
    catch
        eval(['allcells = ',celltype,'cells;']);
    end
    uniquedays = unique(allcells(:,1:2),'rows');  
    
%     load([topDir, 'AnalysesAcrossAnimals\npCells_', region, '.mat'])
%     uniquedays = unique(npCells(:,1:2),'rows');  
%     
    %use value calculated from cs_PDI (from all cells/animals)
    sigPDI = significantPDI.(region)(2);
    
    for ad = 1:size(uniquedays,1)
        
        cells = allcells(ismember(allcells(:,1:2), uniquedays(ad,:), 'rows'),:);
%         cells = npCells(ismember(npCells(:,1:2), uniquedays(ad,:), 'rows'),:);

        if size(cells,1) >= minensemblecells
            
            animal = animals{uniquedays(ad,1)};
            
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            
            day = uniquedays(ad,2);
            daystr = getTwoDigitNumber(day);
            
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            load([animDir, animal,'nosepokeWindow',daystr,'.mat'])
            load([animDir,animal,'spikes',daystr,'.mat'])
            
            runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
            for ep = 1:length(runeps)
                
               
                sessionSpikingData = {};
                epoch = runeps(ep);
                epstr = getTwoDigitNumber(epoch);
                disp(['Doing animal ', animal, ' day ', daystr, ' epoch ', epstr]);
                
                 [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                
                 trigs = odorTriggers{day}{epoch}.allTriggers;
                 npoffs = nosepokeWindow{day}{epoch};
                 
                 switch trialtype
                     case 'correct'
                         lefttrigs = trigs(correct_left);
                         righttrigs = trigs(correct_right);
                     case 'incorrect'
                         lefttrigs = trigs(incorrect_left);
                         righttrigs = trigs(incorrect_right);
                 end
                
                 %trigs are only correct trials, left then right
                trigs = [lefttrigs;righttrigs];
                
                %% Calculate the session PDI using the means
                
                
                for tr = 1:length(trigs)
                    trig = trigs(tr);
                    
                    trigwin = [trig-win(1), trig + win(2)];
                    trigbins = [trig-win(1):binsize:trig + win(2)];
                    
                    %create one large cell array, then use later when
                    %looking at individual trials
                    for c = 1:size(cells,1) 
                            if ~isempty(spikes{day}{epoch}{cells(c,3)}{cells(c,4)}.data)
                                %gather spikes from each trial, separate
                                %into bins
                                cellspikes = spikes{day}{epoch}{cells(c,3)}{cells(c,4)}.data(:,1);
                                winspikes = cellspikes(cellspikes > trigwin(1) & cellspikes <= trigwin(2));
                                binspikecount = histcounts(winspikes',trigbins);

                                sessionSpikingData{c}(tr,:) = binspikecount;
                                
                                
                            else
                                sessionSpikingData{c}(tr,:) = zeros(1, length(bins)-1);
                            end

                    end
                    
                end
                
                %determine how many spikes across all cells occured on each
                %trial. Replace with nans for all trials that did not meet threshold. 
                
                    trialspikes = sum(cell2mat(cellfun(@(x) nansum(x,2), sessionSpikingData,'UniformOutput',false)),2);
                    bad = find(trialspikes < mintrialspikes); 
                for d = 1:length(sessionSpikingData)
                    sessionSpikingData{d}(bad,:) = nan;
                end
                
                
                %convert from spikecount to FR
                sessionSpikingData = cellfun(@(x) x/binsize, sessionSpikingData,'UniformOutput',false);

                %separate into left and right trials
                switch trialtype
                    case 'correct'
                        leftSpikingData = cellfun(@(x) x(1:length(correct_left),:), sessionSpikingData, 'UniformOutput', false);
                        rightSpikingData = cellfun(@(x) x(length(correct_left)+1:end,:), sessionSpikingData, 'UniformOutput', false);
                    case 'incorrect'
                        leftSpikingData = cellfun(@(x) x(1:length(incorrect_left),:), sessionSpikingData, 'UniformOutput', false);
                        rightSpikingData = cellfun(@(x) x(length(incorrect_left)+1:end,:), sessionSpikingData, 'UniformOutput', false);
                end
                
                %remove rows with all nans
                bad = cellfun(@(x) find(all(isnan(x),2)),leftSpikingData,'UniformOutput',false);
                for d = 1:length(leftSpikingData)
                    leftSpikingData{d}(bad{d},:) = [];
                end
                
                bad = cellfun(@(x) find(all(isnan(x),2)),rightSpikingData,'UniformOutput',false);
                for d = 1:length(rightSpikingData)
                    rightSpikingData{d}(bad{d},:) = [];
                end
                
                %get mean PVs across trials
                leftmeanPVs = vertcat(cell2mat(cellfun(@(x) mean(x, 1), leftSpikingData,'UniformOutput',false)'));
                rightmeanPVs = vertcat(cell2mat(cellfun(@(x) mean(x, 1), rightSpikingData,'UniformOutput',false)'));
                
                %smooth mean PVs using interpolation
                leftmeanPVs = smoothdata(interp1(bins(1:end-1),leftmeanPVs',newbins),'gaussian',15)';
                rightmeanPVs = smoothdata(interp1(bins(1:end-1),rightmeanPVs',newbins),'gaussian',15)';

                %% compare individual trials to the session means
                for tr = 1:length(trigs)
                    
                    %get random sample of cells to normalize cell number
                    %across sessions
                    
                    samp  = randsample(size(leftmeanPVs,1),minensemblecells);
                    
                    if tr <= length(lefttrigs)
                        oppositePVs = rightmeanPVs(samp,:);
                    else
                        oppositePVs = leftmeanPVs(samp,:);
                    end
                    
                    trig = trigs(tr);
                    
                        trialspiking = vertcat(cell2mat(cellfun(@(x) x(tr,:), sessionSpikingData(1,samp), 'UniformOutput',false)'));
                        
                        %smooth data using interpolation
                        trialPVs = smoothdata(interp1(bins(1:end-1),trialspiking',newbins),'gaussian',15)';
                        
                        
                        %calc trial PDI
                        trialPDI = zeros(1,length(newbins));
                        for b = 1:length(newbins)
                            trialPV = trialPVs(:,b);
                            oppositePV = oppositePVs(:,b);
                            
                            CCbin = corrcoef(trialPV, oppositePV);
                            
                            PDIbin = 1-abs(CCbin(1,2)); %PV Differential Index, from Igarashi et al 2014
                            if isnan(PDIbin)
                                PDIbin = 0;
                            end
                            
                            trialPDI(b) = PDIbin;
                            
                        end
                        
                        %% determine time at which trial PDI hits the previously-calculated significant PDI value for the whole population
                        %binsaftertrig = find(bins(1:end-1)> 0); %at least 200ms after NP
                        
                        %exclude trial if all bins are significant
                        %moving window of 2 bins, find first time at which
                        %2 consecutive bins are significant (i.e. becomes sig
                        %and stays, not up and down).
                        %if this doesn't happen, exclude trial.

                        timepoint = [];
                        crossings = polyxpoly(newbins, trialPDI, newbins, repmat(sigPDI,1,length(newbins))); 
                        
                        figure
                                plot(newbins, trialPDI,'k-')
                                hold on
                                plot(newbins, repmat(sigPDI,1,length(newbins)), 'r--')
                                plot(crossings,repmat(sigPDI,1,length(crossings)),'gs'); 
                             
                                %find the crossing after which PDI stays
                                %above thresh
                        for t = 1:length(crossings)
                            crossing = crossings(t);
                            prevbin = lookup(crossing,newbins,-1);
                            nextbin = lookup(crossing,newbins,1);
                            
                            %don't use if the previous bin was above
                            %sigPDI (i.e. negative slope), or if it is the last bin, or if it is
                            %within the first 2 bins
                            if trialPDI(prevbin) > sigPDI || nextbin == length(trialPDI) || prevbin == 1 || prevbin == 2
                                continue
                                
                                %must be before the last 2 bins, and all
                                %following PDI must be significant, or if
                                %not the it must be the last positive slope
                                %crossing
                            elseif length(trialPDI) - nextbin <2 || all(trialPDI(nextbin:end-3) > sigPDI) || (t == length(crossings) || t == length(crossings)-1) 
                                
                                timepoint = crossing;
%                                 if timepoint < 0.86 && timepoint > 0.84
%                                     disp('Paused')
%                                 end
                               plot(crossings(t),sigPDI,'m*'); 
                                break
                            end
                        end
%                         
                        
                        if ~isempty(timepoint)
                            trialSigPDI{day}{epoch}(tr,1) = timepoint;
                            trialSigPDI{day}{epoch}(tr,2) = trig;
%                             
                        else
                            trialSigPDI{day}{epoch}(tr,1) = NaN;
                            trialSigPDI{day}{epoch}(tr,2) = NaN;
                        end
                      
                        
                        close all
                end
%                 pause
%                 close all
            end
            %% Save 
                save([animDir,animal,'trialSigPDI_NP_',trialstr,region,daystr,'.mat'],'trialSigPDI')
                clear trialSigPDI
        end
        
    end
    allvals.(region) = vals;
end


 %cs_PDIlfpPowerCorrelation
cs_PDIlfpCoherenceCorrelation_sextiles_NP
cs_PDIBehaviorCorrelation_NP
