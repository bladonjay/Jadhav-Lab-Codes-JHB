[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};
numtrials = [];
for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    epochs = cs_getRunEpochs(animDir, animal, 'odorplace');
    
    days = unique(epochs(:,1));
    for d = 1:length(days)
        day = days(d);
        daystr = getTwoDigitNumber(day);
        
        load([animDir,animal,'odorTriggers',daystr,'.mat']);
        dayeps = epochs(find(epochs(:,1) == day),2); %#ok<FNDSB>
        daytrials = 0;
        for ep = 1:length(dayeps)
            epoch = dayeps(ep);
            eptrials = length(odorTriggers{day}{epoch}.allTriggers);
            daytrials = daytrials + eptrials;
        end
        numtrials = [numtrials; daytrials]; %#ok<AGROW>
       
    end
    
end

histogram(numtrials,[40,50,60,70,80,90,100,110,120,130,140])
xlabel('Number of Trials per Day')
yticks([0 1 2 3 4])
axis([36, 144, 0 3.5])

figtitle = 'Number of Trials Distribution';
figfile = [figDir,'Behavior\',figtitle];
saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);