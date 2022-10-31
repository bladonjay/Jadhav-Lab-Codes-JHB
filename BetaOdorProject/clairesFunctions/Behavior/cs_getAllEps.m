%cs_getAllEps
topDir = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

runtype = 'noodor';
allEps = [];
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    runeps = cs_getRunEpochs(animDir, animal,runtype);
    animvect = repmat(a,size(runeps,1),1);
    allEps = [allEps;animvect,runeps];
end

save([topDir,'AnalysesAcrossAnimals\allEps_',runtype],'allEps');
    
