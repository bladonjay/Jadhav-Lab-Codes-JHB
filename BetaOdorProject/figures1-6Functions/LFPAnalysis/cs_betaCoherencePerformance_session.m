%cs_betaCoherencePerformance

%compare beta power between correct and incorrect trials
%average beta power across nosepoke window, get distribution for each trial
%type

%can also downsample, since there are way fewer incorrect trials
close all
clear
figure, hold on
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%regions = {'CA1-PFC','PFC-OB','CA1-OB'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};

labels_all = [];
results_all = [];
pvals = [];
calc = 1;
figure
for r = 1:length(regions)
    C = [];
    I = [];
    
    region = regions{r};
    disp(['Doing ', region])
    if calc ~=1
        load([topDir,'AnalysesAcrossAnimals\betaCoh_perf'])
        C = betaCoh_perf.(erase(region,'-'))(betaCoh_perf.(erase(region,'-'))(:,2) == 1,1);
        I = betaCoh_perf.(erase(region,'-'))(betaCoh_perf.(erase(region,'-'))(:,2) == 0,1);
    else
        for a = 1:length(animals)
            animal = animals{a};
            disp(['Starting animal ',animal]);
            animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
            daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
            tetinfo = loaddatastruct(animDir, animal,'tetinfo');
            
            days = unique(daymatrix(:,1));
            for day = days'
                daystr = getTwoDigitNumber(day);
                odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
                npWins = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
                
                epochs = daymatrix(daymatrix(:,1) == day,2);
                load([animDir, animal, 'coherence',region,daystr]);
                
                
                for ep = epochs'
                    C_trials = [];
                I_trials = [];
                    epstr = getTwoDigitNumber(ep);
                    
                    time = coherence{day}{ep}.time;
                    %data = coherence{day}{ep}.Coh;
                    data = coherence{day}{ep}.rawCoh;
                    mn = coherence{day}{ep}.rawMean;
                    sd = coherence{day}{ep}.rawSD;
                    
                    %z-score to full epoch
                    data = (data-mn)./sd;
                    
                    freq = coherence{day}{ep}.freq;
                    %find beta freq rows
                    betarows = freq >= 15 & freq <=30;
                    data = data(betarows,:);
                    data = mean(data,1); %mean beta coherence
                    
                    
                    
                    %get trials
                    [cl, cr, il, ir] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    
%                     nptrigs = npWins{day}{ep}(:,2);
%                     nptrigs = [nptrigs-0.5 nptrigs+0.5];
                    nptrigs = npWins{day}{ep};
                    correctwins = nptrigs([cl;cr],:);
                    incorrectwins = nptrigs([il;ir],:);
                    

                    %get beta coherence for each trial
                    for t = 1:size(correctwins)
                        correctbeta = mean(data(isExcluded(time,correctwins(t,:))));
                        C_trials = [C_trials;correctbeta];
                    end
                    for t = 1:size(incorrectwins)
                        incorrectbeta = mean(data(isExcluded(time,incorrectwins(t,:))));
                        I_trials = [I_trials;incorrectbeta];
                    end
                    
                    C_trials = C_trials(~isnan(C_trials));
                I_trials = I_trials(~isnan(I_trials));

                C = [C;mean(C_trials)];
                I = [I;mean(I_trials)];
                end
                
                
            end
            
        end
        %     C = C(~isnan(C));
        %     Cerr = stderr(C);
        %
        %     I = I(~isnan(I));
        %     Ierr = stderr(I);
    end
    
    
    results = [C;I];
    labels = [ones(length(C),1);zeros(length(I),1)];
   betaCoh_perf.(erase(region,'-')) = [results,labels];
    
    p = signrank(C,I);
    disp([region, ' pvalue = ', num2str(p)]);
    
    
    subplot(1,3,r)
    err_C = stderr(C);
    err_I = stderr(I);
    
    cs_errorbar(1,mean(C),err_C,'color','black');
    hold on
    cs_errorbar(2,mean(I),err_I,'color','red');
    
    
    %p = ranksum(C,I);
    pvals = [pvals;p];
    %disp([region, ' pvalue = ', num2str(p)]);
    
    xlim([0 3])
    title(region);
    ylabel('Beta coherence')
    text(0, mean(I),['p = ',num2str(round(p,2,'significant'))]);
    box off
    
    
    
    
end


figfile = [figDir, 'EEG\betaCoh_performance_sessions'];

print('-djpeg', figfile);
print('-dpdf', figfile);

save([topDir,'AnalysesAcrossAnimals\betaCoh_perf_sessions'],'betaCoh_perf');