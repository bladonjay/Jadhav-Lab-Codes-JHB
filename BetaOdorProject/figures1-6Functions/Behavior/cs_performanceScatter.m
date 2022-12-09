%cs_performanceScatter(animals)
[topDir, figDir]= cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

%colors = {'r','g','b','k','m'};
%c = linspace(0,length(animals)+1,length(animals)+1);
sessionmeans = [];
figure, hold on
for a = 1:length(animals)
    animal = animals{a};
    %color = colors{a};
    animperf = [];
    %days = daycell{a};
    PerfAll = [];
    animDir = [topDir,animal,'Expt\', animal, '_direct\'];
    
    runmat = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(runmat(:,1));
    
    for d = 1:length(days)
        day = days(d);
        
        daystring = getTwoDigitNumber(day);
        
        load([animDir,'BinaryPerf\',animal,'BinaryPerf',daystring,'.mat']);
         if strcmpi(animal,'CS41') && day<3
           epochs = 1; % claire collapsed epochs for this animal only
            %dio{1,day}{1}=dio{1,day}; % have to add back in a cell for epoch
        else
            epochs = runmat(runmat(:,1) == day, 2);
        end
        
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            bp = BinaryPerf{day}{epoch};
            perf = sum(bp)/length(bp);
            %mn = mean(bp);
            sessionmeans = [sessionmeans; perf];
            animperf = [animperf; perf];
            %scatter(ones(length(sessionmeans),1), sessionmeans, 'jitter', 'on', 'jitterAmount', 0.1);

        end
    end


%scatter(x,y,25,c(2:6),'filled')
    scatter([ones(length(animperf),1)+(a-1)], animperf, 35, 'filled', 'jitter', 'on', 'jitterAmount', 0.1);

    mn = mean(animperf);
    plot([a-0.15, a+0.15],[mn, mn], 'k-','LineWidth',3);
end

axis([0 length(animals)+1 0 1])
xticklabels(0:length(animals));
ylabel('Fraction Correct')
xlabel('Animal Number')
plot([0 length(animals)+1], [.5 .5], 'k--','LineWidth',2)


figfile = [figDir,'Behavior\BehaviorScatterplot'];
%print('-dpdf', figfile);
%print('-djpeg', figfile);
%print('-dpng', figfile);
%saveas(gcf,figfile,'fig');

