%cs_listPhaseLockedAllCells

%make list of phase locked cells

%% first set up our data
clear
topDir = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
eegregions = {'CA1','PFC','OB'};
regions = {'CA1','PFC'};
trialstr=[];
rhythms={'resp','beta'};
%trialstr='Incorrect';


%% run the code
for rh=1:length(rhythms)
    for cr = 1:length(regions) % cell region
        for er=1:length(eegregions) % eeg region
            
            p1_Allcells=[];
            
            for a = 1:length(animals)
                animal = animals{a};
                animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
                fullDir=[animDir,'PhaseLocking\' trialstr];
                
                fileNames=[animal,'phaselock_', rhythms{rh}, '_', regions{cr}, '-', eegregions{er}, '*'];
                files = dir(fullfile(fullDir,fileNames));
                % there ought to be a single file
                load(fullfile(files(d).folder,files(d).name))
                
                for d = 1:length(length(phaselock))
                    
                    % i already save the variable out as phaselock
                    %eval(['phaselock = ',rhythms{rh},'_phaselock',regions{cr},';']);
                    
                    %                     cellfilter = '($prayl < 0.05) && ($Nspikes > 0)';
                    cellfilter = '$Nspikes > 0';
                    cells = evaluatefilter(phaselock{d},cellfilter);
                    %                     if ~isempty(cells)
                    %                         noeps = cells(:,[2 3]);
                    %                         cells = unique(noeps,'rows');
                    %                         cells = [repmat(a, size(cells,1),1), repmat(day, size(cells,1),1), cells];
                    %
                    %                         plcells = [plcells; cells];
                    %                     end
                    if ~isempty(cells)
                        pl_Allcells = [pl_Allcells; cells];
                    end
                end
            end
            save([dataDir,'PhaseLocking\pl_AllCells_',rhythms{rh},'_',regions{cr},'-',eegregions{er},'.mat'], 'pl_Allcells');
        end
    end
end



