animals = {'CS31','CS33','CS34','CS35'};

regions = 'CA1-PFC';

topDir = 'D:\OdorPlaceAssociation\';

for a = 1:length(animals)
    animal = animals{a};
    
    coherencefiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'coherence',regions,'*']);
    
    for d = 1:length(coherencefiles) 
            load([topDir,animal,'Expt\',animal,'_direct\',coherencefiles(d).name])
            coherence_temp = coherence;
            clear coherence
            
            coherence{1,d} = coherence_temp;
            
            save([topDir,animal,'Expt\',animal,'_direct\',animal, 'coherence',regions, getTwoDigitNumber(d),'.mat'],'coherence');
    end
end

  