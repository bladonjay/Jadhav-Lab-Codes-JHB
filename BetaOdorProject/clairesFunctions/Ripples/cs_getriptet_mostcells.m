
function riptet = cs_getriptet_mostcells(animal,day)

%assigns 'riptet' description to tetinfo for the tetrode that had the most
%cells for each epoch.
topDir = cs_setPaths();
dataDir = [topDir,animal,'Expt',filesep,animal,'_direct',filesep];

load([dataDir,animal,'tetinfo.mat']);
tetfilter = '(isequal($area,''CA1''))';
tets = evaluatefilter(tetinfo{day},tetfilter);
tets = unique(tets(:,2));

epochs = find(~cellfun(@isempty, tetinfo{day}));

numcells = zeros(length(tets),1);
for ep = 1:length(epochs)
    epoch = epochs(ep);

    for t = 1:length(tets)
        tet = tets(t);
        if isfield(tetinfo{day}{epoch}{tet}, 'descrip')
            if strcmp(tetinfo{day}{epoch}{tet}.descrip, 'riptet')
                tetinfo{day}{epoch}{tet} = rmfield(tetinfo{day}{epoch}{tet},'descrip');
            end
        end
        
        if ~isempty(tetinfo{day}{epoch}{tet})
            num = tetinfo{day}{epoch}{tet}.numcells;
        else
            num = 0;
        end
        numcells(t) = numcells(t)+num;
    end
    
    
end
[~,ind] = max(numcells);
riptet = tets(ind);