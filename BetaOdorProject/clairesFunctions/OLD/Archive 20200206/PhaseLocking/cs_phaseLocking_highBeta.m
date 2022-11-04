%calculates phase locking for each epoch, stores as beta_phaselocking file
%in animal's _direct folder, rather than combined in AnalysesAcrossAnimals.
%
%uses high beta times, calculated from cs_getHighBetaTimes and stored as
%"highBeta" in animal_dir folder.
clear
[topDir, figDir] = cs_setPaths;
gethighbeta = 1;
if gethighbeta == 0
    cs_getHighBetaTimes;
end

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

regions = {'CA1','PFC','OB'};
freq = 'beta';
trialtypes = {'incorrect','correct'};
%trialtypes = {'prelearn','postlearn'};

for tt = 1:length(trialtypes)
    trialtype = trialtypes{tt};
    %trialtype = 'correct';
    
    switch trialtype
        case 'correct'
            trialstr = '';
        case 'incorrect'
            trialstr = '_incorrect';
        case 'prelearn'
            trialstr = '_prelearn';
        case 'postlearn'
            trialstr = '_postlearn';  
    end
    
    for r = 1:length(regions)
        region = regions{r};
        for a = 1:length(animals)
            animal = animals{a};
            
            animDir = ([topDir, animal,'Expt\',animal, '_direct\']);

            daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');

            load([animDir,animal,'highBeta',trialstr,'.mat']);
            %         load([animDir,animal,'betaWindows.mat']);
            %         highBeta = betaWindows;
            %trigfiles = dir();
            days = unique(daymatrix(:,1));
            
            for d = 1:length(days)
                day = days(d);
                daystr = getTwoDigitNumber(day);
                disp(['Doing ', animal,' Day ', daystr]);
                
                if strcmp(trialstr,'_incorrect')
                    cs_calcPhaseLocking_beta_incorrect(a, day, region, highBeta)
                else
                    cs_calcPhaseLocking_beta(a, day, region, highBeta, trialstr)
                end 
                %                 case 'theta'
                        %                     cs_calcPhaseLocking_beta(a, day, region, highBeta, trialstr)
                
                
            end
        end
        cs_listPhaseLockedCells(freq,region,trialstr)
    end
end