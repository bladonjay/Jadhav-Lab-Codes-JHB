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
calc = 0;
bootstrap = 1;
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
                epstr = getTwoDigitNumber(ep);
                
                time = coherence{day}{ep}.time;
                data = coherence{day}{ep}.rawCoh;
                mn = coherence{day}{ep}.rawMean;
                sd = coherence{day}{ep}.rawSD;
                freq = coherence{day}{ep}.freq;
                
                %z-score to full epoch
                data = (data-mn)./sd;
                
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
        p =   min([p_upper, p_lower]);      
        pvals = [pvals;p];
        
        subplot(1,3,r)
        
        CI = getCI(Cdist,95);
        plot([1 1], CI','k-')
        hold on
        plot(1,mean(Cdist),'ko')
        plot(2,mean(I),'ro');

        xlim([0 3])
        
            switch region
                case 'CA1-PFC'
                    ylim([0.24 0.32])
                case 'CA1-OB'
                    ylim([0.18 0.27])
                case 'PFC-OB'
                    ylim([-0.09 0.01])
            end
        %ylim([-0.1 0.3])
        title(region);
        ylabel('beta coherence')
        text(0, mean(I),['p = ',num2str(round(p,2,'significant'))]);
        
        results = [C;I];
        labels = [ones(length(C),1);zeros(length(I),1)];
        betaCoh_perf.(erase(region,'-')) = [results,labels];
        box off
                disp([region, ' pvalue = ', num2str(p)]);
%          if length(region) <7
%             regionstr = [region,' '];
%         else
%             regionstr = region;
%          end
%         
%         errorbar(r, mean(C), stderr(C),'b.')
%         hold on
%         errorbar(r+0.25, mean(Idist), 2*std(Idist),'r.')
%          
%         labels = [repmat([regionstr,'c'],length(C),1); repmat([regionstr,'i'],length(Idist),1)];
%         results_all = [results_all;results];
%         labels_all = [labels_all;labels];
%         
%         if length(results_all) ~= length(labels_all)
%             keyboard
%         end
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

figfile = [figDir, 'EEG\betaCoh_performance'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile); 

    save([topDir,'AnalysesAcrossAnimals\betaCoh_perf'],'betaCoh_perf');
% xticks([1.125 2.125 3.125]);
% xticklabels(regions);
% 
% axis([0.75 3.75 -0.2 0.4]);
% text(1.5, -.2, ['p = ', num2str(round(pvals(1),2))]);
% text(2.5, -.2, ['p = ', num2str(round(pvals(2),2))]);
% text(3.5, -.2, ['p = ', num2str(round(pvals(3),2))]);
% ylabel('Beta (15-30 Hz) Power');
% figfile = [figDir, 'EEG\BetaCoherence_performance'];
%             
%     print('-djpeg', figfile);
%     print('-dpdf', figfile);
% 
% 
% figure
% hold on
% boxplot(results_all,labels_all,'Colors','br','symbol','k.');
% xticks([1.5 3.5 5.5 7.5]);
% xticklabels(regions);
% text(1.5, 4, ['p = ', num2str(round(pvals(1),2))]);
% text(3.5, 4, ['p = ', num2str(round(pvals(2),2))]);
% text(5.5, 4, ['p = ', num2str(round(pvals(3),2))]);
% 
% 
% % text(1,0.16,['p = ', num2str(p)]);
% %plot([0 r+1.25],[0 0],'k--');
% ylabel('Z-scored Beta Coherence');
% % xticks(1:length(regions))
% % xticklabels(regions);
% legend({'Correct', 'Incorrect'},'Location','Northeast')
% 
% figfile = [figDir, 'Behavior\betaCoherence_performance'];
%             
%     print('-djpeg', figfile);
%     print('-dpdf', figfile);
