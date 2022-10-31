%cs_betaZscore

%gets beta power, zscored to baseline period- 1s before odor onset.

close all
clear
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};

baseline = 'session'; %session or pre 

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    nosepokeWindows = loaddatastruct(animDir, animal,'nosepokeWindow');
    tetinfo = loaddatastruct(animDir, animal,'tetinfo');
    dayepochs = cs_getRunEpochs(animDir, animal,'odorplace');
    days = unique(dayepochs(:,1));
    
    disp(['Doing animal ',animal]);
    for r = 1:length(regions)
        region = regions{r};
        
        for day = days'
            epochs = dayepochs(dayepochs(:,1) == day,2);
            
            for ep = epochs'
                %get tetrode with most cells
                tetfilter = ['isequal($area,''',region,''')'];
                tets = evaluatefilter(tetinfo{day}{ep},tetfilter);
                numcells = cellfetch(tetinfo{day}{ep}(tets),'numcells');
                [~,ind] = max(cell2mat(numcells.values));
                tet = tets(ind(1));
                tetstr = getTwoDigitNumber(tet);
                beta = loadeegstruct(animDir, animal, 'beta',day,ep,tet);
                betatime = geteegtimes(beta{day}{ep}{tet});
                beta = double(beta{day}{ep}{tet}.data(:,3));
                
                %get baseline
                switch baseline
                    case 'session'
                        baselineMean = mean(beta);
                        baselineStd = std(beta);
                    case 'pre'
                        prestimtimes = [nosepokeWindows{day}{ep}(:,1)-1, nosepokeWindows{day}{ep}(:,1)];
                        prestimbeta = beta(isExcluded(betatime,prestimtimes));
                        baselineMean = mean(prestimbeta);
                        baselineStd = std(prestimbeta);
                end
                
                %get zscore for each trial
                stimtimes = nosepokeWindows{day}{ep};
                alltrigbeta = [];
                for t = 1:size(stimtimes,1)
                    trig = stimtimes(t,:);
                    trigbeta = beta(isExcluded(betatime,trig));
                    alltrigbeta = [alltrigbeta;mean(trigbeta)];
                end
                
                zscorebeta = (alltrigbeta-baselineMean)/baselineStd;
                betaZscore{day}{ep}.(region) = zscorebeta;
            end
            
        end
        
    end
    save([animDir,animal,'betaZscore'],'betaZscore');
    clear betaZscore
    
end