function cs_addriptet(animal)

%assigns 'riptet' description to tetinfo for the tetrode that had the most
%cells for each epoch. 
topDir = cs_setPaths();
dataDir = [topDir,animal,'Expt',filesep,animal,'_direct',filesep];

load([dataDir,animal,'tetinfo.mat']);
tetfilter = '(isequal($area,''CA1''))';

days = find(~cellfun(@isempty, tetinfo));

for d = 1:length(days)
    day = days(d);
    epochs = find(~cellfun(@isempty,tetinfo{day}));
    
    for ep = 1:length(epochs)
        epoch = epochs(ep);
        
        tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
        
        numcells = [];
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
            numcells = [numcells;num];
        end
        [~,ind] = max(numcells);
        riptet = tets(ind);
        tetinfo{day}{epoch}{riptet}.descrip = 'riptet';
        
    end
end
save([dataDir,animal,'tetinfo.mat'],'tetinfo');