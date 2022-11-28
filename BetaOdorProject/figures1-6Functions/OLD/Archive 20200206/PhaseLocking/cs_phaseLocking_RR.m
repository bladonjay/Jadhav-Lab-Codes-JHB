%calculates phase locking for each epoch, stores as beta_phaselocking file
%in animal's _direct folder, rather than combined in AnalysesAcrossAnimals.
%
%uses high beta times, calculated from cs_getHighBetaTimes and stored as
%"highBeta" in animal_dir folder.

[topDir, figDir] = cs_setPaths;

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
freq = 'resp';
trialtype = 'incorrect'; 

switch trialtype
    case 'correct'
        trialstr = '';
    case 'incorrect'
        trialstr = '_incorrect';
end

for r = 1:length(regions)
    region = regions{r};
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        epochs = cs_getRunEpochs(animDir, animal, 'odorplace');
        days = unique(epochs(:,1));
        
        %just take entire nosepoke window since the RR amplitude is fairly
        %consistent throughout odor period
        npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
        odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            disp(['Doing ', animal,' Day ', daystr]);
            if strcmp(trialstr,'_incorrect')
                    cs_calcPhaseLocking_RR_incorrect(a, day, region, npWins)
                    
            else
                    cs_calcPhaseLocking_RR(a, day, region, npWins, trialstr)
            end 
                
            
           
        end
    end
    cs_listPhaseLockedCells(freq,region,trialstr)
end
