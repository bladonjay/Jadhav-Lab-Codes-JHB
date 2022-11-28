%Calculates and saves average firing rate during odor on each trial for each
%odor-responsive cell

clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
topDir = cs_setPaths();

    %get nosepoke cells
    load([topDir,'AnalysesAcrossAnimals\npCells_CA1']);
    allCells= npCells;
    load([topDir, 'AnalysesAcrossAnimals\npCells_PFC']);
    allCells= [allCells;npCells];


for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        spikes = loaddatastruct(animDir, animal,'spikes');
        nosepokeWindow = loaddatastruct(animDir, animal,'nosepokeWindow');
        
        cells = allCells(allCells(:,1) == a,(2:4));
        
        days = unique(cells(:,1));
        
        for d = 1:length(days)
            day = days(d);

            daycells = cells(cells(:,1) == day,:);
            
            runEps = cs_getRunEpochs(animDir, animal,'odorplace',day);
            epochs = unique(runEps(:,2));
            for c = 1:size(daycells,1)
                cell = daycells(c,:);

                rates = [];
                for ep = epochs'
%                     diffs = diff(nosepokeWindow{day}{ep}');
%                     half = diffs'./2;
%                     wins = [nosepokeWindow{day}{ep}(:,2)-half,nosepokeWindow{day}{ep}(:,2)];
                    wins = nosepokeWindow{day}{ep};
                    
                     spiketimes = spikes{cell(1)}{ep}{cell(2)}{cell(3)}.data;
                    
                     %get firing rate on each trial
                     for w = 1:size(wins)
                         
                     trigspikes = spiketimes(isExcluded(spiketimes, wins(w,:)));
                     trialfr = length(trigspikes)/(wins(w,2)-wins(w,1));
                     
                     rates = [rates;trialfr];
                     end
                     
                     %collect into cell array
                    nosepokeFR{day}{ep}{cell(2)}{cell(3)} = rates;
                end
            end
        end
        
        %save in animal directory
        save([animDir, animal,'nosepokeFR'],'nosepokeFR');
        clear nosepokeFR
end
                
              