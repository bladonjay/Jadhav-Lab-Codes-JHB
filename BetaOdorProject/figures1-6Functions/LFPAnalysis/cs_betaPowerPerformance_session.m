%cs_betaPowerPerformance

%compare beta power between correct and incorrect trials

%get means across sessions (days) and compare these distributions directly
%rather than downsampling/bootstrapping


close all
clear
[topDir,figDir] = cs_setPaths();
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};
figure
hold on
labels_all = [];
results_all = [];
pvals = [];
calc = 0;
for r = 1:length(regions)
    C = [];
    I = [];
    
    region = regions{r};
    disp(['Doing ', region])
    if calc ~=1
        load([topDir,'AnalysesAcrossAnimals\betaPower_perf_sessions'])
        C = betaPower_perf.(region)(betaPower_perf.(region)(:,2) == 1,1);
        I = betaPower_perf.(region)(betaPower_perf.(region)(:,2) == 0,1);
    else
        for a = 1:length(animals)
            animal = animals{a};
            disp(['Doing animal ',animal]);
            animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
            daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
            tetinfo = loaddatastruct(animDir, animal,'tetinfo');
            odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers');
            npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
            rewards = loaddatastruct(animDir, animal,'rewards');
            days = unique(daymatrix(:,1));
            for day = days'
                daystr = getTwoDigitNumber(day);
                epochs = daymatrix(daymatrix(:,1) == day,2);
                
                %             tetfilt = ['strcmp($area,''',region,''')'];
                %             tets = evaluatefilter(tetinfo,tetfilt);
                
                %load([animDir, animal, 'specfloor',region,daystr]);
                C_trials = [];
                I_trials = [];
                for ep = epochs'
                    epstr = getTwoDigitNumber(ep);
                    
                    tet = cs_getMostCellsTet(animal, day, ep, region);
                    if isempty(tet)
                        break
                    end
                    
                    beta = loadeegstruct(animDir, animal, 'beta',day,ep,tet);
                    time = geteegtimes(beta{day}{ep}{tet});
                    data = double(beta{day}{ep}{tet}.data(:,3));
                    
                    %get baseline for zscore
                    rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
                    %rtimes = sort(rtimes(:,1));
                    
                    rewarddata = mean(data(isExcluded(time,rtimes)));
                    %mn = mean(rewarddata);
                    sd = std(data(isExcluded(time,rtimes)));;
                    
                    %                 mn = mean(data);
                    %                 sd = std(data);
                    %zscore
                    data = (data-rewarddata)/sd;
                    
                    [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    correctwins = npWins{day}{ep}([correct_left;correct_right],:);
                    incorrectwins = npWins{day}{ep}([incorrect_left;incorrect_right],:);
                    
                    %get beta power for each trial
                    for t = 1:size(correctwins)
                        correctbeta = mean(data(isExcluded(time,correctwins(t,:))));
                        C_trials = [C_trials;correctbeta];
                    end
                    for t = 1:size(incorrectwins)
                        incorrectbeta = mean(data(isExcluded(time,incorrectwins(t,:))));
                        I_trials = [I_trials;incorrectbeta];
                    end
                end
                
                C = [C;mean(C_trials)];
                I = [I;mean(I_trials)];
                
            end
            
        end
%         C = C(~isnan(C));
%         Cerr = stderr(C);
%         
%         I = I(~isnan(I));
%         Ierr = stderr(I);
        
    end
    
    
    % switch region
    %     case 'CA1'
    %         ylim([-0.05 0.05])
    %     case 'PFC'
    %         ylim([0.20 0.33])
    %     case 'OB'
    %         ylim([0.86 1.15])
    % end
    
    results = [C;I];
    labels = [ones(length(C),1);zeros(length(I),1)];
    betaPower_perf.(region) = [results,labels];
    
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
    ylabel('Beta power')
    text(0, mean(I),['p = ',num2str(round(p,2,'significant'))]);
    box off
    
    
    
    
end


figfile = [figDir, 'EEG\BetaPower_performance_sessions'];

print('-djpeg', figfile);
print('-dpdf', figfile);

save([topDir,'AnalysesAcrossAnimals\betaPower_perf_sessions'],'betaPower_perf');
% % figure
% % hold on
% % boxplot(results_all,labels_all,'Colors','br','symbol','k.');
% % xticks([1.5 3.5 5.5 7.5]);
% % xticklabels(regions);
% % text(1.5, 4, ['p = ', num2str(round(pvals(1),2))]);
% % text(3.5, 4, ['p = ', num2str(round(pvals(2),2))]);
% % text(5.5, 4, ['p = ', num2str(round(pvals(3),2))]);
% % %text(7.5, 7, ['p = ', num2str(round(pvals(4),2))]);
% % %axis([0 r+1.25 -0.1 1.1])
% % %text(1,0.16,['p = ', num2str(p)]);
% % %plot([0 r+1.25],[0 0],'k--');
% % ylabel('Beta (15-30 Hz) Power');
% % %legend({'Correct', 'Incorrect'},'Location','East')
% %
% % figfile = [figDir, 'EEG\BetaPower_performance_boxplot'];
% %
% %     print('-djpeg', figfile);
% %     print('-dpdf', figfile);
