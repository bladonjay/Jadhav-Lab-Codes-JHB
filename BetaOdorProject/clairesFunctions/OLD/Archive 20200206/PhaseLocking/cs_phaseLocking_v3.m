%calculates phase locking for each epoch, stores as beta_phaselocking file
%in animal's _direct folder, rather than combined in AnalysesAcrossAnimals.
%

[topDir, figDir] = cs_setPaths;

animals = {'CS31','CS33','CS34','CS35'};

for a = 1:length(animals)
    animal = animals{a};
    
    animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
    
    trigfiles = dir([animDir,animal,'odorTriggers*']);
    numdays = length(trigfiles);
    
     for d = 1:numdays
        load(trigfiles(d).name);
        day = length(odorTriggers);
        daystr = getTwoDigitNumber(day);
        disp(['Doing ', animal,' Day ', daystr]);
            
        cswt_beta_phaselocking(animal, day, odorTriggers)
        
     end
end

cs_listPhaseLockedCells