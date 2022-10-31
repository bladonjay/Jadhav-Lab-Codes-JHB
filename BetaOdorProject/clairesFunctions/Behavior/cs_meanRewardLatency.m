clear
animals = {'CS31','CS33','CS34','CS35','CS44'};
topDir = cs_setPaths;

latencies = [];
for a = 1:length(animals)
    animal = animals{a};

animDir = [topDir,animal,'Expt\',animal,'_direct\'];
rewards = loaddatastruct(animDir, animal,'rewards');
odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');

runeps = cs_getRunEpochs(animDir, animal,'odorplace');
days = unique(runeps(:,1));
    
rtimes = [];
for day = days'
    eps = runeps((runeps(:,1) == day),2);
    for ep = eps'
        lreward = rewards{day}{ep}.leftWindows(:,1);
        rreward = rewards{day}{ep}.rightWindows(:,1);
        allrew = sort([lreward;rreward]);
        
        [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
        
        lodor = odorTriggers{day}{ep}.allTriggers(cl);
        rodor = odorTriggers{day}{ep}.allTriggers(cr);
        odors = sort([lodor;rodor]);
        r = [lreward-lodor; rreward-rodor];
        rtimes = [rtimes;r];
    end
end
meantime = mean(rtimes);
latencies = [latencies;meantime];
end

latency = mean(latencies);
save([topDir, 'AnalysesAcrossAnimals\rewardLatency'],'latency');
