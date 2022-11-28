%compare average beta power and coherence between odor disengagement and
%reward disengagement

clear
close all
regions = {'CA1-PFC'};
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
win = [0.5 0.5];

[topDir,figDir] = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];

freqs = {'beta'};

for f = 1:length(freqs)
    freq = freqs{f};
    
    figure
    for r = 1:length(regions)
        region = regions{r};
        
        allsessions = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir,animal,'Expt\', animal, '_direct\'];
            
            daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
            days = unique(daymatrix(:,1));
            
            %load data
            odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
            npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
            rewards = loaddatastruct(animDir, animal,'rewards');
            
            for day = days'
                daystr = getTwoDigitNumber(day);
                epochs = daymatrix(daymatrix(:,1) == day,2);
                load([animDir, animal, 'coherence', region, daystr]);

                for ep = epochs'
                    epstr = getTwoDigitNumber(ep);
                    
                    cohfull = coherence{day}{ep}.Coh;
                    time = coherence{day}{ep}.time;
                    switch freq
                        case 'beta'
                            freqrange = [15 30];
                        case 'resp'
                            freqrange = [7 8];
                    end
                    rows = find(coherence{day}{ep}.freq >= freqrange(1) & coherence{day}{ep}.freq <= freqrange(2));
                    data = mean(cohfull(rows,:),1);
                    
%                         epmean = mean(data);
%                         epsd = std(data);
%                         
%                         data = (data-epmean)/epsd;
                    [cl, cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    nptrigs = sort(npWins{day}{ep}([cl;cr],2));
                    nptrigs = [nptrigs-win(1), nptrigs+win(2)];
                    rewardtrigs = sort([rewards{day}{ep}.leftWindows(:,2); rewards{day}{ep}.rightWindows(:,2)]);
                    rewardtrigs = [rewardtrigs-win(1), rewardtrigs+win(2)];
                    
                    npdata = data(isExcluded(time,nptrigs));
                    rewdata = data(isExcluded(time,rewardtrigs));
                    
                    allsessions = [allsessions; mean(npdata), mean(rewdata)];
                    
                end
            end
        end
        
        p = signrank(allsessions(:,1),allsessions(:,2));
        np_mn = mean(allsessions(:,1));
        np_err = stderr(allsessions(:,1));
        rew_mn = mean(allsessions(:,2));
        rew_err = stderr(allsessions(:,2));
        
        subplot(1,length(regions),r)
        cs_errorbar(1,np_mn,np_err,'color','black');
        %plot([1,1], [np_mn-np_err, np_mn+np_err], 'k-');
        hold on
        cs_errorbar(2,rew_mn,rew_err,'color','blue');
        %plot([1.25,1.25], [rew_mn-rew_err, rew_mn+rew_err], 'b-');
        text(1,rew_mn+rew_err,['p = ',num2str(p)])
        xlabel(region)
        xlim([0 3])
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        ylabel(['Z-scored ', freq, ' coherence']);
        box off
       
    end
%     axP = get(gca,'Position');
%     fP = get(gcf,'Position');
    legend({'odor exit','','reward exit'},'Location','NorthEastOutside')
%     set(gca,'Position',axP);
%     fP(3) = fP(3)+200;
%     set(gcf,'Position',fP);

    figfile = [figDir,'EEG\','coherenceOdorVsReward_',freq];
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
end