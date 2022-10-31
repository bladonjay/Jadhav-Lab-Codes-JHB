%cs_avgBetaOnset
clear
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
avgBetaOnset = [];
for a = 1:length(animals)
    animal = animals{a};
    animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
    
    betaWindows = loaddatastruct(animDir,animal,'betaWindows');
    days = 1:length(betaWindows);
    for day = days
        npWindows = loaddatastruct(animDir,animal,'nosepokeWindow',day);
        
        epochs = cs_getRunEpochs(animDir, animal, 'odorplace',day);
        for ep = epochs(:,2)'
            
            wins = npWindows{day}{ep}(:,1);
            betastart = betaWindows{day}{ep}(:,1);
            
            starttimes = wins - betastart;
            
            avgstart = nanmean(starttimes);
            
            avgBetaOnset = [avgBetaOnset;avgstart];
        end
        
    end
    
end

avgBetaOnset = nanmean(avgBetaOnset);