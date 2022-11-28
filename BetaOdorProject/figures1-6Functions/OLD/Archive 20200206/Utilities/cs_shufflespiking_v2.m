function [leftfr, rightfr] = cs_shufflespiking_v2(lefttrials, righttrials, cellinds,binsize) 

%does one iteration, gives new mean FR

%load in spiking data for each cell on each trial, along with cellindicies
cellinds = cell2mat(cellinds);
animals = unique(cellinds(:,1));
leftfr = []; rightfr = [];
for a = 1:length(animals)
    animal = animals(a);
    
    animcells = cellinds((cellinds(:,1) == animal),2:end);
    
    days = unique(animcells(:,1));
    
    for d = 1:length(days)
        day = days(d);
        daycells = animcells((animcells(:,1) == day),2:end);
        daycells = [repmat(animal,size(daycells,1),1), repmat(day,size(daycells,1),1), daycells];
        
        [~,daycellinds] = ismember(daycells, cellinds, 'rows');

        numlefttrials = size(lefttrials{daycellinds(1)},1);
        numrighttrials = size(righttrials{daycellinds(1)},1);
        
        shuffledindex = randperm(numlefttrials+numrighttrials)';
        
        for c = 1:size(daycells,1)
            left = lefttrials{daycellinds(c)};
            right = righttrials{daycellinds(c)};
            alltrials = [left; right];
            shuffledtrials = alltrials(shuffledindex,:);
            
            newleft = shuffledtrials(1:numlefttrials,:);
            newright = shuffledtrials(numlefttrials+1:end,:);
            
            leftfr = [leftfr; mean(newleft,1)/binsize];
            rightfr = [rightfr; mean(newright,1)/binsize];
        end
    end
end