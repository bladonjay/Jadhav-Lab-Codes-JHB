animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

[topDir, figDir]= cs_setPaths();

%colors = {'r','g','b','k','m'};
%c = linspace(0,length(animals)+1,length(animals)+1);
sessionmeans_o = [];
sessionmeans_a = [];
%figure, hold on
 odorperf = [];
 airperf = [];
for a = 1:length(animals)
    animal = animals{a};
    %color = colors{a};
    %days = daycell{a};
    
    animDir = [topDir,animal,'Expt\', animal, '_direct\'];
    
    runmat_odor = cs_getRunEpochs(animDir, animal, 'odorplace');
    runmat_air = cs_getRunEpochs(animDir,animal,'noodor');
    
    if isempty(runmat_air)
        continue
    end
    
    days = unique(runmat_odor(:,1));
    
    for d = 1:length(days)
        day = days(d);
        
        daystring = getTwoDigitNumber(day);
        
        load([animDir,'BinaryPerf\',animal,'BinaryPerf',daystring,'.mat']);
        
        epochs = runmat_odor(runmat_odor(:,1) == day, 2);
        session = [];
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            bp = BinaryPerf{day}{epoch};
            perf = sum(bp)/length(bp);
            %mn = mean(bp);
            session = [session; perf];
            odorperf = [odorperf; perf];
            %scatter(ones(length(sessionmeans),1), sessionmeans, 'jitter', 'on', 'jitterAmount', 0.1);

        end
        sessionmeans_o = [sessionmeans_o;mean(session)];
    end
    
    days = unique(runmat_air(:,1));
    for d = 1:length(days)
        day = days(d);
        
        daystring = getTwoDigitNumber(day);
        
        load([animDir,'BinaryPerf\',animal,'BinaryPerf',daystring,'.mat']);
        
        epochs = runmat_air(runmat_air(:,1) == day, 2);
        session = [];
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            bp = BinaryPerf{day}{epoch};
            perf = sum(bp)/length(bp);
            %mn = mean(bp);
            session = [session; perf];
            airperf = [airperf; perf];
            %scatter(ones(length(sessionmeans),1), sessionmeans, 'jitter', 'on', 'jitterAmount', 0.1);

        end
         sessionmeans_a = [sessionmeans_a;mean(session)];
    end
    
end

figure
%     scatter(ones(1,length(odorperf)), odorperf, 35, 'filled', 'jitter', 'on', 'jitterAmount', 0.1);
%     hold on
%     scatter(ones(1,length(airperf))+1,airperf,35,'filled','jitter', 'on', 'jitterAmount', 0.05);
% xlim([0 3])
%     mn = mean(odorperf);
%     plot([0.85, 1.15],[mn, mn], 'k-','LineWidth',3);
% 
%     mn = mean(airperf);
%     plot([1.85 2.15],[mn, mn], 'k-','LineWidth',3);
%     
%     p = ranksum(odorperf,airperf);

scatter(ones(1,length(sessionmeans_o)), sessionmeans_o, 35, 'filled', 'jitter', 'on', 'jitterAmount', 0.1);
    hold on
    scatter(ones(1,length(sessionmeans_a))+1,sessionmeans_a,35,'filled','jitter', 'on', 'jitterAmount', 0.05);
xlim([0 3])
    mn = mean(sessionmeans_o);
    plot([0.85, 1.15],[mn, mn], 'k-','LineWidth',3);

    mn = mean(sessionmeans_a);
    plot([1.85 2.15],[mn, mn], 'k-','LineWidth',3);
    
    p = ranksum(sessionmeans_o,sessionmeans_a);
    
    
ylabel('Fraction Correct')
xticks([1 2]);
xticklabels({'Odor Sessions','Air Sessions'})
text(1,0.5,['p = ',num2str(p)])
ylim([0 1]);

figfile = [figDir,'Behavior\OdorVsAirPerf'];
print('-dpdf', figfile);
print('-djpeg', figfile);
%print('-dpng', figfile);
%saveas(gcf,figfile,'fig');
