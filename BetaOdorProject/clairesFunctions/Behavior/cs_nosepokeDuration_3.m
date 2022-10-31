%cs_nosepokeDuration_3

%finds average np duration across SESSIONS

[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

duration_all = [];
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    
    nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow');
    odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    
    for day = days'
        epochs = runeps(runeps(:,1)==day,2);
        
        for ep = epochs'
            
            %correct trials
            [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
            nps = nosepokeWindow{day}{ep};
            duration = nps(:,2)-nps(:,1);
            
            duration_all = [duration_all;mean(duration)];
            
            
        end
    end
end
mn = mean(duration_all);
se = stderr(duration_all);
%sd = std(duration_all);
