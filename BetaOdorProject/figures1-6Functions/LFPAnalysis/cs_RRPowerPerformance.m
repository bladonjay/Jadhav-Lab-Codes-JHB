%cs_betaPowerPerformance

%compare beta power between correct and incorrect trials
%average beta power across nosepoke window, get distribution for each trial
%type

%can also downsample, since there are way fewer incorrect trials
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
bootstrap = 1;
for r = 1:length(regions)
    C = [];
    I = [];
    
    region = regions{r};
    disp(['Doing ', region])
    if calc ~=1
        load([topDir,'AnalysesAcrossAnimals\RRPower_perf'])
        C = RRPower_perf.(region)(RRPower_perf.(region)(:,2) == 1,1);
        I = RRPower_perf.(region)(RRPower_perf.(region)(:,2) == 0,1);
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
                for ep = epochs'
                    epstr = getTwoDigitNumber(ep);
                    
                    tet = cs_getMostCellsTet(animal, day, ep, region);
                    if isempty(tet)
                        break
                    end
                    
                    resp = loadeegstruct(animDir, animal, 'resp',day,ep,tet);
                    time = geteegtimes(resp{day}{ep}{tet});
                    data = double(resp{day}{ep}{tet}.data(:,3));
                    
                    
                    %get baseline for zscore - use reward times
                rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
               
            
                rewarddata = data(isExcluded(time,rtimes));
                mn = mean(rewarddata);
                sd = std(rewarddata);
                
                
%                     mn = mean(data);
%                     sd = std(data);
                    %zscore
                    data = (data-mn)/sd;
                    
                    [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                    correctwins = npWins{day}{ep}([correct_left;correct_right],:);
                    incorrectwins = npWins{day}{ep}([incorrect_left;incorrect_right],:);
                    
                    %get beta power for each trial
                    for t = 1:size(correctwins)
                        correctresp = mean(data(isExcluded(time,correctwins(t,:))));
                        C = [C;correctresp];
                    end
                    for t = 1:size(incorrectwins)
                        incorrectresp = mean(data(isExcluded(time,incorrectwins(t,:))));
                        I = [I;incorrectresp];
                    end
                end
                
            end
            
        end
        C = C(~isnan(C));
        Cerr = stderr(C);
        
        I = I(~isnan(I));
        Ierr = stderr(I);
        
    end
    
    if bootstrap == 1
        iterations = 1000;
        Cdist = [];
        for i = 1:iterations
            samp = datasample(1:size(C,1),size(I,1));
            trials = C(samp);
            Cdist = [Cdist;mean(trials)];
        end
        p_upper  = sum(Cdist >= mean(I)) /length(Cdist);
        p_lower = sum(Cdist <= mean(I))/length(Cdist);
        p = min([p_upper, p_lower]);
        pvals = [pvals;p];
        
        subplot(1,3,r)
        %get CI
        CI = getCI(Cdist,95);
        plot([1 1], CI','k-')
        hold on
        plot(1,mean(Cdist),'ko')
        plot(2,mean(I),'ro');
        
        xlim([0 3])
        title(region);
        ylabel('RR power')
        text(0, mean(I),['p = ',num2str(round(p,2,'significant'))]);
        
        switch region
            case 'CA1'
                ylim([0.25 0.45])
            case 'PFC'
                ylim([0.75 1.0])
            case 'OB'
                ylim([2.0 2.5])
        end
        
        results = [C;I];
        labels = [ones(length(C),1);zeros(length(I),1)];
        RRPower_perf.(region) = [results,labels];
        box off
        disp([region, ' pvalue = ', num2str(p)]);
        %text(r-0.1,0.2,['p = ', num2str(p)]);
        
        %         errorbar(r, mean(C), stderr(C),'b.')
        %         hold on
        %         errorbar(r+0.25, mean(Cdist), 2*std(Cdist),'r.')
        %
        %plot(r, mean(C),'k.');
        %axis([0 r+1.25 0.1 0.75])
        %         data{r}.C = mean(C);
        %         data{r}.I = Idist;
        
        %         results = [I;Cdist];
        %         labels = [repmat([region(1),'c'],length(C),1); repmat([region(1),'i'],length(Cdist),1)];
        %         results_all = [results_all;results];
        %         labels_all = [labels_all;labels];
        %         %plotvars = [plotvars; results,labels];
        
    else
        %calc significance
        %[~,p] = ttest2(C,I);
        p = ranksum(C,I);
        pvals = [pvals;p];
        disp([region, ' pvalue = ', num2str(p)]);
        
        %plot
        %         errorbar(r,mean(C),Cerr,'k.');
        %         hold on
        %         errorbar(r+0.25,mean(I),Ierr,'k.');
        %         data{r}.C = C;
        %         data{r}.I = I;
        results = [C;I];
        labels = [repmat([region(1),'c'],length(C),1); repmat([region(1),'i'],length(I),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        
    end
    
    
    %rrPowerData.(region) = C;
    %save([topDir,'AnalysesAcrossAnimals\rrPowerData'],'rrPowerData');
end

figfile = [figDir, 'EEG\RRPower_performance'];

print('-djpeg', figfile);
print('-dpdf', figfile);

save([topDir,'AnalysesAcrossAnimals\RRPower_perf'],'RRPower_perf');

% xticks([1.125 2.125 3.125 4.125]);
% xticklabels(regions);
%
% axis([0.75 3.75 -0.2 1.2]);
% text(1.5, -.2, ['p = ', num2str(round(pvals(1),2))]);
% text(2.5, -.2, ['p = ', num2str(round(pvals(2),2))]);
% text(3.5, -.2, ['p = ', num2str(round(pvals(3),2))]);
% %text(4.5, -.2, ['p = ', num2str(round(pvals(4),2))]);

% ylabel('Respiratory Rhythm (7-8 Hz) Power');
% figfile = [figDir, 'EEG\RespPower_performance'];
%
%     print('-djpeg', figfile);
%     print('-dpdf', figfile);
%
% figure
% hold on
% boxplot(results_all,labels_all,'Colors','br','symbol','k.');
% xticks([1.5 3.5 5.5 7.5]);
% xticklabels(regions);
% text(1.5, 4, ['p = ', num2str(round(pvals(1),2))]);
% text(3.5, 4, ['p = ', num2str(round(pvals(2),2))]);
% text(5.5, 4, ['p = ', num2str(round(pvals(3),2))]);
% %text(7.5, 7, ['p = ', num2str(round(pvals(4),2))]);
% %axis([0 r+1.25 -0.1 1.1])
% %text(1,0.16,['p = ', num2str(p)]);
% %plot([0 r+1.25],[0 0],'k--');
% ylabel('Respiratory Rhythm (7-8 Hz) Power');
% %legend({'Correct', 'Incorrect'},'Location','East')
%
% figfile = [figDir, 'EEG\RespPower_performance_boxplot'];
%
%     print('-djpeg', figfile);
%     print('-dpdf', figfile);
