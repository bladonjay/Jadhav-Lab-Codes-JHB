%How many selective cells are place cells? 
%Overlap
function [cells, fraction, num] = cs_odorSelectivePlaceCells()

topDir = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
%figDir = [figDir,'NicePPTFigures\'];
%freqbands = {'beta','theta'};
regions = {'CA1'};

for r = 1:length(regions)
    region = regions{r};
    load([dataDir, 'selectiveCells_',region,'.mat']);
    load([dataDir, 'placeCells_',region,'.mat']);

    cells = intersect(selectivecells,placecells,'rows');
    
    fraction = size(cells,1) / size(selectivecells,1);
    num = size(cells,1);
    
end

end
