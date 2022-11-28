%cs_getOdorOffset
topDir = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal, '_direct\'];
    npWins = loaddatastruct(animDir, animal,'nosepokeWindow');
    
    days = find(~cellfun(@isempty,npWins));
    for d = days
        daystr = getTwoDigitNumber(d);
        eps = find(~cellfun(@isempty,npWins{d}));
        
        for ep = eps
            odorOffset{d}{ep} = npWins{d}{ep}(:,2);
        end
        
        save([animDir,animal,'odorOffset',daystr],'odorOffset');
        clear odorOffset
    end

end