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
bootstrap = 1;
for r = 1:length(regions)
    C = [];
    I = [];

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
            odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
            npWins = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
            
            epochs = daymatrix(daymatrix(:,1) == day,2);
            load([animDir, animal, 'coherence',region,daystr]);
            
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                time = coherence{day}{ep}.time;
                data = coherence{day}{ep}.Coh;
                freq = coherence{day}{ep}.freq;
                
                %find beta freq rows
                betarows = freq >= 15 & freq <=30;
                data = data(betarows,:);
                data = mean(data,1); %mean beta coherence
                
                
                %get trials
                [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                correctwins = npWins{day}{ep}([correct_left;correct_right],:);
                incorrectwins = npWins{day}{ep}([incorrect_left;incorrect_right],:);
                
                %get beta coherence for each trial
                for t = 1:size(correctwins)
                    correctbeta = mean(data(isExcluded(time,correctwins(t,:))));
                    C = [C;correctbeta];
                end
                for t = 1:size(incorrectwins)
                    incorrectbeta = mean(data(isExcluded(time,incorrectwins(t,:))));
                    I = [I;incorrectbeta];
                end
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
        pvals = [pvals;p];
        
        results = [C;Idist];
        
         if length(region) <7
            regionstr = [region,' '];
        else
            regionstr = region;
         end
        
        errorbar(r, mean(C), stderr(C),'b.')
        hold on
        errorbar(r+0.25, mean(Idist), 2*std(Idist),'r.')
         
        labels = [repmat([regionstr,'c'],length(C),1); repmat([regionstr,'i'],length(Idist),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        
        if length(results_all) ~= length(labels_all)
            keyboard
        end
        %plotvars = [plotvars; results,labels];
    
    else 
         %calc significance
        [~,p] = ttest2(C,I);
        p = ranksum(C,I);
        pvals = [pvals;p];
        disp([region, ' pvalue = ', num2str(p)]);
        
        if length(region) <7
            regionstr = [region,' '];
        else
            regionstr = region;
        end
       
        results = [C;I];
        labels = [repmat([regionstr,'c'],length(C),1); repmat([regionstr,'i'],length(I),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        
    end
    
end

xticks([1.125 2.125 3.125]);
xticklabels(regions);

axis([0.75 3.75 -0.2 0.4]);
text(1.5, -.2, ['p = ', num2str(round(pvals(1),2))]);
text(2.5, -.2, ['p = ', num2str(round(pvals(2),2))]);
text(3.5, -.2, ['p = ', num2str(round(pvals(3),2))]);
ylabel('Beta (15-30 Hz) Power');
figfile = [figDir, 'EEG\BetaCoherence_performance'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);


figure
hold on
boxplot(results_all,labels_all,'Colors','br','symbol','k.');
xticks([1.5 3.5 5.5 7.5]);
xticklabels(regions);
text(1.5, 4, ['p = ', num2str(round(pvals(1),2))]);
text(3.5, 4, ['p = ', num2str(round(pvals(2),2))]);
text(5.5, 4, ['p = ', num2str(round(pvals(3),2))]);


% text(1,0.16,['p = ', num2str(p)]);
%plot([0 r+1.25],[0 0],'k--');
ylabel('Z-scored Beta Coherence');
% xticks(1:length(regions))
% xticklabels(regions);
legend({'Correct', 'Incorrect'},'Location','Northeast')

figfile = [figDir, 'Behavior\betaCoherence_performance'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);
