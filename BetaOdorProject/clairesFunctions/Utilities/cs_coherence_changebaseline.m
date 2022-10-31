%coherence_changebaseline

topDir = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};

for r = 1:length(regions)
    region = regions{r};
    
    for a = 1:length(animals)
        animal = animals{a};
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        rewards = loaddatastruct(animDir, animal,'rewards');
        daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
        
        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = daymatrix(daymatrix(:,1) == day,2);
            load([animDir, animal, 'coherence',region,daystr]);
            for ep = epochs'
                time = coherence{day}{ep}.time;
                data = coherence{day}{ep}.Coh;
                freq = coherence{day}{ep}.freq;
                sd = coherence{day}{ep}.sd;
                mn = coherence{day}{ep}.mean;
                
                %-- un-zscore
                rawCoh = (data.*sd)+mn;
                
                %-- get baseline based on reward times
                 rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
               
            
                rewarddata = rawCoh(isExcluded(time,rtimes));
                mn = mean(rewarddata);
                sd = std(rewarddata);
                
                %-- new zscore
                
                data = (rawCoh-mn)./sd;
                
                out.Coh = data;
                out.rawCoh = rawCoh;
                out.rawMean = mean(rawCoh,2);
                out.rawSD = std(rawCoh,0,2);
                out.time = time;
                out.freq = freq;
                out.mean_base = mn;
                out.sd_base = sd;
               
                
                coherence{day}{ep} = out;
                
            end
            save([animDir,animal,'coherence',region,daystr],'coherence')
            clear coherence
        end
    end
end