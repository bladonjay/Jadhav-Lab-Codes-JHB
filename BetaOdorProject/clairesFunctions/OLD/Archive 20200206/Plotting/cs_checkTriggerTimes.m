function cs_checkTriggerTimes(dataDir, animal)

for d = 1:length(dir([animal,'odorTriggers*']))
    
    daystring = getTwoDigitNumber(d);
    
    load([dataDir,animal,'odorTriggers',daystring,'.mat'])

    for i = 1:length(odorTriggers{1,d})
        if ~isempty(odorTriggers{1,d}{1,i})
            all = odorTriggers{1, d}{1, i}.allTriggers;
            correct = odorTriggers{1, d}{1, i}.correctTriggers;
            incorrect = odorTriggers{1, d}{1, i}.incorrectTriggers;

            figure, hold on
            plot(all, 1, 'ko')
            plot(correct, 1, 'gx')
            plot(incorrect, 1, 'rx')
            
        end
    end
end

