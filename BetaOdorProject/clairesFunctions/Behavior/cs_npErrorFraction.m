%cs_npErrorFraction
animals = {'CS31','CS33','CS34','CS35','CS39','CS42','CS44'};
topDir = cs_setPaths();
allerrfract = [];
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    runEps = cs_getRunEpochs(animDir, animal,'odorplace');
    days = unique(runEps(:,1));
    npErrors = loaddatastruct(animDir, animal,'npErrors');
    nosepokeWindow = loaddatastruct(animDir, animal,'nosepokeWindow');
    
    for day = days'
        eps = runEps(runEps(:,1) ==day,2);
        errs = 0;
        good = 0;
        for ep = eps'
            errs = errs+ length(npErrors{day}{ep});
            good = good + size(nosepokeWindow{day}{ep},1);
        end
        
        total = errs +good;
        
        fract = errs/total;
        allerrfract = [allerrfract;fract];
    end
        
    
end

mn = mean(allerrfract);
sem = stderr(allerrfract);

success = 1-mn;
