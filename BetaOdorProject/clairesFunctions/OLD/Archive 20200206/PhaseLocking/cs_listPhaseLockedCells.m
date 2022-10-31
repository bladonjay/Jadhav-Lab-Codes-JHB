%make list of phase locked cells
function cs_listPhaseLockedCells(freq,eegregion,trialstr)
topDir = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%betaregion = 'OB';
regions = {'CA1','PFC'};


for r = 1:length(regions)
    region = regions{r};
    plcells = [];
    
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        files = dir([animDir,'PhaseLocking\',animal,freq,'phaselock_',region,'-',eegregion,trialstr,'_0*']);
        
        for d = 1:length(files)
            load([animDir, 'PhaseLocking\',files(d).name])
            eval(['phaselock = ',freq,'_phaselock',region,';']);
            day = length(phaselock);
            
            cellfilter = '($prayl < 0.05) && ($Nspikes > 0)';
            cells = evaluatefilter(phaselock{day},cellfilter);
       
            if ~isempty(cells)
                noeps = cells(:,[2 3]);
                cells = unique(noeps,'rows');
                cells = [repmat(a, size(cells,1),1), repmat(day, size(cells,1),1), cells];

                plcells = [plcells; cells];
            end
        end
    end
    
    save([dataDir,'PhaseLocking\plCells_',freq,'_',region,'-',eegregion,trialstr,'.mat'], 'plcells');
end


        