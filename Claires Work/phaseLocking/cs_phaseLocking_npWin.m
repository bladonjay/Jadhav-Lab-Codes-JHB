%calculates phase locking for each epoch, stores as beta_phaselocking file
%in animal's _direct folder, rather than combined in AnalysesAcrossAnimals.
%

[topDir, figDir] = cs_setPaths;

animals = {'CS31','CS33','CS34','CS35'};
eegregions = {'CA1','PFC','OB'};
freq = 'theta';

for r = 1:length(eegregions)
    eegregion = eegregions{r};
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        
        trigfiles = dir([animDir,animal,'nosepokeWindow*']);
        numdays = length(trigfiles);
        
        for d = 1:numdays
            load(trigfiles(d).name);
            day = length(nosepokeWindow);
            daystr = getTwoDigitNumber(day);
            disp(['Doing ', animal,' Day ', daystr]);
            
            switch freq
                case 'beta'
                    cs_calcPhaseLocking_beta(animal, day, eegregion, nosepokeWindow)
                case 'theta'
                    cs_calcPhaseLocking_theta(animal, day, eegregion, nosepokeWindow)
            end
            
        end
    end
    cs_listPhaseLockedCells(freq,eegregion)
end
