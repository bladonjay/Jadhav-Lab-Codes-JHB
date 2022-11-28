% Claire Symanski - March 2 2017
% Takes Binary performance files and combines them for specified days, and
% outputs statespace plot

%cs_stateSpacePerf_plotAllAnimals({'CS31','CS33','CS34','CS35'}, {[1,2,3],[1,2,3,4],[1,2,3,4],[1,2,3]})


function cs_performanceBar(animals)
[topDir, figDir]= cs_setPaths();
%daycell should be cell array with matrix of days to use for each animal.
%I.e. daycell = {[1,2,3],[1,2,3,4],[1,2,3,4]};
%in same order as animals ({'CS31','CS33','CS34'})
%cs_stateSpacePerf_plotAllAnimals({'CS31','CS33','CS34','JS11'}, {[1,2,3],[1,2,3,4],[1,2,3,4],[1,2,3,4,5]})

%
%colors = {[rgb('HotPink')], [rgb('MediumSlateBlue')], [rgb('MediumTurquoise')], [rgb('Tomato')], [rgb('LimeGreen')], [rgb('HotPink')]};

meanall = [];
for a = 1:length(animals)
    animal = animals{a};
    
    
    %days = daycell{a};
    PerfAll = [];
    animDir = [topDir,animal,'Expt\', animal, '_direct\'];
    
    runmat = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(runmat(:,1));
    
    for d = 1:length(days)
        day = days(d);
        
        daystring = getTwoDigitNumber(day);
        
        load([animDir,'BinaryPerf\',animal,'BinaryPerf',daystring,'.mat']);
        
        epochs = runmat(runmat(:,1) == day, 2);
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            bp = BinaryPerf{day}{epoch};
            
            PerfAll = [PerfAll; bp];
        end
    end


meanperf= mean(PerfAll);
meanall(a) = meanperf;
stderror(a) = std(PerfAll) / sqrt( length(PerfAll));
stdev(a) = std(PerfAll);
end

figure,
xvals = 1:length(animals);
bar(xvals, meanall,'FaceColor',rgb('Black'));
hold on
% h(1) = bar(1, meanall(1),'FaceColor',colors{1}); hold on
% h(2) = bar(2, meanall(2),'FaceColor',colors{2});
% h(3) = bar(3, meanall(3),'FaceColor',colors{3});
% h(4) = bar(4, meanall(4),'FaceColor',colors{4});
% h(5) = bar(5, meanall(5),'FaceColor',colors{5});

plot([0 length(animals)+1], [.5 .5], 'r--','LineWidth',3)

errorbar(xvals,meanall,stderror,'.','LineWidth',3,'Color',rgb('DimGrey'))

%legend('Rat 1','Rat 2','Rat 3','Rat 4','Rat 5', 'Location', 'eastoutside')
box off

axis([0 length(animals)+1 0 1])
%xticks([])
yticks([0 .25 .5 .75 1])

set(gcf, 'Position', [2000 250 600 500]);
set(gca,'fontsize',20);
box off
%xlabel('Animal Number','FontSize',20)
set(gca,'XColor',[0 0 0]);


%ylabel({'Mean Performance'}, 'FontSize',20)
set(gca,'YColor',[0 0 0]);

%figdir = 'E:\Figures\NicePPTFigures\';
figdir = [figDir,'Behavior\'];
%cd(figdir)

%daystring = num2str(days);

figfile = [figdir,'BehaviorBar'];
print('-dpdf', figfile);
print('-djpeg', figfile);
%print('-dpng', figfile);
saveas(gcf,figfile,'fig');

end