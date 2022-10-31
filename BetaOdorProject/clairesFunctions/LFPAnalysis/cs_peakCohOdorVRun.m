close all
clear
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS44'};
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
        
        [~,peakodorind] = max(odor,[],2);
        [~,peakrunind] = max(run,[],2);
        
        peakodor = newfreq(peakodorind);
        peakrun = newfreq(peakrunind);

        p_CA1odorRun = signrank(peakodor,peakrun);
        
        cs_boxplot([peakodor';peakrun'],[ones(length(peakodor),1);ones(length(peakrun),1)+1])

        %ylim([4 11])
        ylabel('Peak Coherence Frequency (Hz)')
        xticklabels({'Odor','Run'})

        
        text(7,0.1,['p = ',num2str(p_CA1odorRun)])


        figfile = [figDir,'EEG\CohOdorVrun_boxplot'];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        
        
    end


