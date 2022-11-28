clear
close all

[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS33'};
regions = {'CA1','PFC'};
lfpregions = {'CA1','PFC','OB'};
freqs = {'beta','resp'};

for f = 1:length(freqs)
    freq = freqs{f};
    
    for L = 1:length(lfpregions)
        lfpregion = lfpregions{L};
        
        for r = 1:length(regions)
            region = regions{r};
            
            allLFPPower = [];
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
                        
                        for e = 1:length(epochs)
                            epoch = epochs(e);
                            
                            [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                            [correctinds,order] = sort([correct_left;correct_right]);
                            windows = nosepokeWindow{day}{epoch}(correctinds,:);
                            
                            %get neural divergence times, re-sort so it matches
                            %with correct trial order
                            divTimes = trialSigPDI{day}{epoch}(order);
                            
                            %get lfp Power
                            lfptet = cs_getMostCellsTet(animal,day,epoch,lfpregion);
                            
                            lfp = loadeegstruct(animDir, animal, freq,day,epoch,lfptet);
                            lfp = lfp{day}{epoch}{lfptet};
                            times = geteegtimes(lfp);
                            lfp = double(lfp.data(:,3));
                            
                            lfpbins = periodAssign(times, windows); %Assign spikes to align with each trials(same number = same trial, number indicates trial)
                            goodeeg = lfp(find(lfpbins));
                            lfpbins = nonzeros(lfpbins);
                            bp = [];
                            for s = unique(lfpbins)'
                                binpower = mean(goodeeg(lfpbins == s));
                                bp(s,1) = binpower;
                            end
                            
                            %zscore power
                            zscoredpower = (bp-mean(lfp))/std(lfp);
                            
                            if length(bp)<length(divTimes)
                                test = unique(lfpbins,'stable');
                                %find time windows that were not found in lfp
                                notfound = setxor(test,[1:length(divTimes)]);
                                divTimes(notfound) = [];
                            end
                            allLFPPower = [allLFPPower; zscoredpower];
                            allDivTimes = [allDivTimes; divTimes];
                            
                            
                            %allPeakTimes = [allPeakTimes; peakTimes];
                        end
                    end
                end
            end
            
            
            keepinds = ~isnan(allDivTimes);
            allDivTimes = allDivTimes(keepinds);
            allLFPPower = allLFPPower(keepinds);
            
            %% correlation
            
            figure
            plot(allDivTimes, allLFPPower, 'k.','MarkerSize',12)
            hold on
            
            
            fit = polyfit(allDivTimes, allLFPPower,1);
            plot([min(allDivTimes), max(allDivTimes)], polyval(fit,[min(allDivTimes), max(allDivTimes)] ))
            
            %set(gcf,'Position',[1200,300,510,420]);
            xlabel('First significant divergence time')
            ylabel([freq,' power'])
            title([region,'PV - ',lfpregion,' ',freq])
            [CCdiv,p_div] = corrcoef(allDivTimes,allLFPPower);
            R = CCdiv(1,2);
            p = p_div(1,2);
            
            legend(['R = ',num2str(R)], ['p = ', num2str(p)]);
            %text(0.2,1,['R = ',num2str(R), newline, 'p = ', num2str(p)])
            
            figtitle = ['PVDivergence',region,'-',lfpregion,freq,'Power'];
            figfile = [figDir,'PopSelectivity\',figtitle];
            %saveas(gcf,figfile,'fig');
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            
        end
    end
end
