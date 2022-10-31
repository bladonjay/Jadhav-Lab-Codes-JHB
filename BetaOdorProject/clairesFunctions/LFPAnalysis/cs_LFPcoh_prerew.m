close all
clear
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};

labels_all = [];
results_all = [];
pvals = [];
calc = 0;
freqs = {'beta'};
for f = 1:length(freqs)
    fr = freqs{f};
    for r = 1:length(regions)
        pre = [];
        rew = [];
        
        region = regions{r};
        disp(['Doing ', region])
        
        for a = 1:length(animals)
            animal = animals{a};
            disp(['Starting animal ',animal]);
            animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
            
            daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
            tetinfo = loaddatastruct(animDir, animal,'tetinfo');
            odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
            npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
            rewards = loaddatastruct(animDir, animal, 'rewards');
            
            days = unique(daymatrix(:,1));
            for day = days'
                daystr = getTwoDigitNumber(day);
                
                
                epochs = daymatrix(daymatrix(:,1) == day,2);
                load([animDir, animal, 'coherence',region,daystr]);
                eppre = [];
                eprew = [];
                for ep = epochs'
                    epstr = getTwoDigitNumber(ep);
                    
                    time = coherence{day}{ep}.time;
                    data = coherence{day}{ep}.rawCoh;
                    
                    freq = coherence{day}{ep}.freq;
                    
                    %z-score to full epoch
                    
                    rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
                    switch fr
                        case 'beta'
                            %find beta freq rows
                            rows = freq >= 15 & freq <=30;
                            mn = coherence{day}{ep}.rawMean;
                            sd = coherence{day}{ep}.rawSD;
                            data = (data-mn)./sd;
                        case 'resp'
                            rows = freq >= 7 & freq <=8;
                            
                            rewarddata = data(isExcluded(time,rtimes));
                            mn = mean(rewarddata);
                            sd = std(rewarddata);
                            data = (data-mn)./sd;
%                             mn = coherence{day}{ep}.rawMean;
%                             sd = coherence{day}{ep}.rawSD;
%                             data = (data-mn)./sd;
                    end
                    
                    
                    data = data(rows,:);
                    data = mean(data,1); %mean beta coherence
                    
                    [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    wins = npWins{day}{ep}([cl;cr],:);
                    dur = wins(:,2)-wins(:,1);
                    prewins = [wins(:,1)-dur, wins(:,1)];
                    rewwins = [rtimes(:,1),rtimes(:,1)+dur];
                    %get beta coherence for each trial
                    for t = 1:size(prewins,1)
                        
                        predata = mean(data(isExcluded(time,prewins(t,:))));
                        rewdata = mean(data(isExcluded(time,rewwins(t,:))));
                        
                        if ~isnan(predata) && ~isnan(rewdata)
                            
                            eppre = [eppre;predata];
                            eprew = [eprew;rewdata];
                        end
                    end
                    
                end
%                 pre = [pre;mean(eppre)];
%                 post = [post;mean(eppost)];

                     pre = [pre;eppre];
                rew = [rew;eprew];
                
            end
            
        end
        
        
        
        mnpre = mean(pre);
        mnrew = mean(rew);
        errpre = stderr(pre);
        errrew = stderr(rew);
        p = signrank(pre,rew);
        
        figure
        errorbar([1,2],[mnpre, mnrew],[errpre,errrew])
        text(1,mnpre,['p = ',num2str(p)]);
        xticks([1 2])
        xticklabels({'pre','reward'});
        xlim([0 3])
        
        figfile = [figDir,'EEG\CohPreReward_',fr,'_',region];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        
        
    end
end


