%move phaselocking files to new folder, "PhaseLocking" in each animal's
%directory

[topDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    cd(animDir);
    
    files = dir([animal, '*phaselock*']);
    
    if ~exist([animDir,'PhaseLocking'], 'dir')
       mkdir([animDir,'PhaseLocking'])
    end
    cd([animDir,'PhaseLocking'])
    for f = 1:length(files)
        file = [animDir,files(f).name];
        movefile(file);
    end
end

        
    
    
    