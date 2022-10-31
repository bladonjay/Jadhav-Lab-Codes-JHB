clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
regions = {'CA1','PFC'};
lfpregions = {'CA1-PFC'};
freqs = {'beta','resp'};

for f = 1:length(freqs)
    freq = freqs{f};
    switch freq
    case 'beta'
        bandpass = [15 30];
    case 'resp'
        bandpass = [7 9];
    end
    
    
    for L = 1:length(lfpregions)
        lfpregion = lfpregions{L};
        figure
        rnum = 0;
            for r = 1:length(regions)
            region = regions{r};
            rnum = rnum+1;
            allLFPCoh = [];
            allDivTimes = [];
            for a = 1:length(animals)
                animal = animals{a};
                
                animDir = [topDir,animal,'Expt\',animal,'_direct\'];
                
                files = dir([animDir,animal,'trialSigPDI_',region,'*']);
                nosepokeWindow = loaddatastruct(animDir,animal,'nosepokeWindow');
                odorTriggers = loaddatastruct(animDir,animal,'odorTriggers');
                if ~isempty(files)
                    for f = 1:length(files)
                        load([animDir,files(f).name]);
                        day = length(trialSigPDI);
                        daystr = getTwoDigitNumber(day);
                        epochs = find(~cellfun(@isempty,trialSigPDI{day}));
                        
                        load([animDir, animal,'coherence',lfpregion,daystr]);
                        
                        for e = 1:length(epochs)
                            epoch = epochs(e);
                            
                            [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                            [correctinds,order] = sort([correct_left;correct_right]);
                            windows = nosepokeWindow{day}{epoch}(correctinds,:);
                            
                            %get neural divergence times, re-sort so it matches
                            %with correct trial order
                            divTimes = trialSigPDI{day}{epoch}(order);
                            
                            
                            lfp = coherence{day}{epoch}.Coh;
                            goodrows = coherence{day}{epoch}.freq >= bandpass(1) & coherence{day}{epoch}.freq <= bandpass(2);
                            lfp = lfp(find(goodrows),:);
                            lfp = mean(lfp,1);
                            
                            times = coherence{day}{epoch}.time;
                            
                            lfpbins = periodAssign(times, windows); %Assign spikes to align with each trials(same number = same trial, number indicates trial)
                            goodlfp = lfp(find(lfpbins));
                            lfpbins = nonzeros(lfpbins);
                            bp = [];
                            for s = unique(lfpbins)'
                                binpower = mean(goodlfp(lfpbins == s));
                                bp(s,1) = binpower;
                            end
                            
                            %zscore power
                            %zscoredpower = (bp-mean(lfp))/std(lfp);
                            
                            if length(bp)<length(divTimes)
                                test = unique(lfpbins,'stable');
                                %find time windows that were not found in lfp
                                notfound = setxor(test,[1:length(divTimes)]);
                                divTimes(notfound) = [];
                            end
                            allLFPCoh = [allLFPCoh; bp];
                            allDivTimes = [allDivTimes; divTimes];
                            
                            
                            %allPeakTimes = [allPeakTimes; peakTimes];
                        end
                    end
                end
            end
            
            
            keepinds = ~isnan(allDivTimes);
            allDivTimes = allDivTimes(keepinds);
            allLFPCoh = allLFPCoh(keepinds);
            
            %% correlation
            
            subplot(1,2,rnum)
            plot(allDivTimes, allLFPCoh, 'k.','MarkerSize',12)
            hold on
            
            
            fit = polyfit(allDivTimes, allLFPCoh,1);
            plot([min(allDivTimes), max(allDivTimes)], polyval(fit,[min(allDivTimes), max(allDivTimes)] ))
            
            %set(gcf,'Position',[1200,300,510,420]);
            xlabel('PV divergence time')
            ylabel([freq,' coherence'])
            title([region,'PV - ',lfpregion,' ',freq])
            [CCdiv,p_div] = corrcoef(allDivTimes,allLFPCoh);
            R = CCdiv(1,2);
            p = p_div(1,2);
            
            legend(['R = ',num2str(R)], ['p = ', num2str(p)]);
            %text(0.2,1,['R = ',num2str(R), newline, 'p = ', num2str(p)])
            title([region,'PV-',lfpregion,])
            
            
        end
    end
    figtitle = ['PVDivergence_',freq,'Coherence'];
            figfile = [figDir,'PopSelectivity\',figtitle];
            %saveas(gcf,figfile,'fig');
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            
end
