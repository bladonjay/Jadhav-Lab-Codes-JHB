%cs_betaPowerPerformance

%compare beta power between correct and incorrect trials
%uses zscored beta power compared to baseline (1s before odor)

%can also downsample, since there are way fewer incorrect trials
close all
clear
figure
hold on
[topDir,figDir] = cs_setPaths();
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};

labels_all = [];
results_all = [];
pvals = [];
bootstrap = 0;
for r = 1:length(regions)
    C = [];
    I = [];

    region = regions{r};
    disp(['Doing ', region])
    for a = 1:length(animals)
        animal = animals{a};
        disp(['Doing animal ',animal]);
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
        tetinfo = loaddatastruct(animDir, animal,'tetinfo');
        odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers');

        
        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            betaZscore = loaddatastruct(animDir, animal, 'betaZscore');
            
            epochs = daymatrix(daymatrix(:,1) == day,2);
            
            for ep = epochs'
                
                beta = betaZscore{day}{ep}.(region);
                [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                correctbeta = beta([correct_left;correct_right]);
                incorrectbeta = beta([incorrect_left;incorrect_right]);
                
                C = [C;correctbeta];
                I = [I;incorrectbeta];

            end
            
        end
        
    end
    C = C(~isnan(C));
    Cerr = stderr(C);
    
    I = I(~isnan(I));
    Ierr = stderr(I);
    
     if bootstrap == 1
        iterations = 1000;
        Idist = [];
        for i = 1:iterations
            samp = datasample(1:size(I,1),size(C,1));
            trials = I(samp);
            Idist = [Idist;mean(trials)];
        end
        p = sum(Idist >= mean(C)) /length(Idist);
       
        results = [C;Idist];
        labels = [repmat([region(1),'c'],length(C),1); repmat([region(1),'i'],length(Idist),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        %plotvars = [plotvars; results,labels];
    
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

    betaPowerData.(region) = C;
    save([topDir,'AnalysesAcrossAnimals\betaPowerData_zscore'],'betaPowerData');
end
figure
hold on
boxplot(results_all,labels_all,'Colors','br','symbol','w.');
xticks([1.5 3.5 5.5 7.5]);
xticklabels(regions);
text(1.5, 4, ['p = ', num2str(round(pvals(1),2))]);
text(3.5, 4, ['p = ', num2str(round(pvals(2),2))]);
text(5.5, 4, ['p = ', num2str(round(pvals(3),2))]);

axis([0 7 -2 7])
%text(1,0.16,['p = ', num2str(p)]);
%plot([0 r+1.25],[0 0],'k--');
ylabel('Z-scored Beta Power');
%xticks(1:length(regions))
%xticklabels(regions);
legend({'Correct', 'Incorrect'},'Location','East')

figfile = [figDir, 'Behavior\betaPower_performance_zscore'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);
