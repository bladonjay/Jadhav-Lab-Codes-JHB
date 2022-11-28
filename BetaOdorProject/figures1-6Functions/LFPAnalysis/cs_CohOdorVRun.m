close all
clear
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1-PFC'};

labels_all = [];
results_all = [];
pvals = [];
%calc = 0;
%freqs = {'resp'};
    for r = 1:length(regions)
        odor = [];
        run = [];
        
        region = regions{r};
        disp(['Doing ', region])
        
        for a = 1:length(animals)
            animal = animals{a};
            disp(['Starting animal ',animal]);
            animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
            
            daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
            tetinfo = loaddatastruct(animDir, animal,'tetinfo');
            
            days = unique(daymatrix(:,1));
            for day = days'
                daystr = getTwoDigitNumber(day);
                npWins = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
                rewards = loaddatastruct(animDir, animal, 'rewards',day);
                
                epochs = daymatrix(daymatrix(:,1) == day,2);
                load([animDir, animal, 'coherence',region,daystr]);
                epodor = [];
                eprun = [];
                for ep = epochs'
                    epstr = getTwoDigitNumber(ep);
                    
                    time = coherence{day}{ep}.time;
                    data = coherence{day}{ep}.rawCoh;
                    
                    freq = coherence{day}{ep}.freq;
                    stepsize = freq(2)-freq(1);
                    
                    %z-score to full epoch
                    rows = freq >= 6 & freq <=12;
                            mn = coherence{day}{ep}.rawMean;
                            sd = coherence{day}{ep}.rawSD;
                            data = (data-mn)./sd;
%                        
%                     
%                     switch fr
%                         case 'beta'
%                             %find beta freq rows
%                             rows = freq >= 15 & freq <=30;
%                             mn = coherence{day}{ep}.rawMean;
%                             sd = coherence{day}{ep}.rawSD;
%                             data = (data-mn)./sd;
%                         case 'resp'
%                             rows = freq >= 7 & freq <=8;
%                             rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
%                             rewarddata = data(isExcluded(time,rtimes));
%                             mn = mean(rewarddata);
%                             sd = std(rewarddata);
%                             data = (data-mn)./sd;
% %                             mn = coherence{day}{ep}.rawMean;
% %                             sd = coherence{day}{ep}.rawSD;
% %                             data = (data-mn)./sd;
%                     end
                    freq = freq(rows);
                    newfreq = freq(1):stepsize/5:freq(end);
                    
                    data = data(rows,:);
                    %data = mean(data,1); %mean beta coherence
                    
                    wins = npWins{day}{ep};
                    dur = wins(:,2)-wins(:,1);
                    runwins = [wins(:,2), wins(:,2)+dur];
                    %get beta coherence for each trial
                    for t = 1:size(wins,1)
                        
                        rundata = mean(data(:,isExcluded(time,runwins(t,:))),2);
                        odordata = mean(data(:,isExcluded(time,wins(t,:))),2);
                        
                        if ~any(isnan(rundata)) && ~any(isnan(odordata))
                            
                            odordata = smoothdata(interp1(freq,odordata,newfreq),'gaussian',15);
                            rundata = smoothdata(interp1(freq,rundata,newfreq),'gaussian',15);
                            
                            epodor = [epodor;odordata];
                            eprun = [eprun;rundata];
                        end
                    end
                    
                end
                odor = [odor;mean(epodor)];
                run = [run;mean(eprun)];

%                      pre = [pre;eppre];
%                 post = [post;eppost];
                
            end
            
        end
        
        mnodor = mean(odor);
        mnrun = mean(run);
        
        norm_odor = rescale(mnodor);
        norm_run = rescale(mnrun);
        
        figure
        hold on
        plot(newfreq, norm_odor);
        plot(newfreq, norm_run);
        
       [~,peakodorind] = max(odor,[],2);
       [~,peakrunind] = max(run,[],2);
       
        peakodor = newfreq(peakodorind);
        peakrun = newfreq(peakrunind);
        
        p = signrank(peakodor,peakrun);
        
        [maxodor,ind] = max(norm_odor);
plot([newfreq(ind) newfreq(ind)],[0 maxodor],'k--');

[maxrun,ind] = max(norm_run);
plot([newfreq(ind) newfreq(ind)],[0 maxrun],'k--');

text(7,0.1,['p = ',num2str(p)])
ylabel('Normalized CA1-PFC Coherence')
xlabel('Frequency (Hz)');

        figfile = [figDir,'EEG\CohOdorVrun'];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        
        
    end


