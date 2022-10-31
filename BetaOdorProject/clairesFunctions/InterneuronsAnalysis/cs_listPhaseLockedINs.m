%make list of phase locked cells

topDir = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\PhaseLocking\'];
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
eegregions = {'CA1','PFC','OB'};
freqs = {'beta','resp'};


for f = 1:length(freqs)
    freq = freqs{f};
    for er = 1:length(eegregions)
        eegregion = eegregions{er};
        
        for cr = 1:length(cellregions)
            cellregion = cellregions{cr};
            
            plcells = [];
            
            for a = 1:length(animals)
                animal = animals{a};
                animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
                
                try
                load([animDir,'PhaseLocking\Interneurons\',animal,'phaselock_',freq,'_',cellregion,'-',eegregion,'.mat'])
                catch
                    continue
                end
                
                cellfilter = '($prayl < 0.05) && ($Nspikes > 0)';
                cells = evaluatefilter(phaselock,cellfilter);
                
                if ~isempty(cells)
                    noeps = cells(:,[1 3 4]);
                    cells = unique(noeps,'rows');
                    cells = [repmat(a, size(cells,1),1),cells];
                    
                    plcells = [plcells; cells];
                end
            end
            
            
            save([dataDir,'plINs_',freq,'_',cellregion,'-',eegregion,'.mat'], 'plcells');
        end
    end
end