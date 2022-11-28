function tet = cs_getMostCellsTet(animal,day,epoch,region)

topDir = cs_setPaths();
animDir = [topDir, animal,'Expt\',animal,'_direct\'];
tetinfo = loaddatastruct(animDir, animal, 'tetinfo');

filt = ['strcmp($area,''',region,''')'];
tets = evaluatefilter(tetinfo{day}{epoch},filt);
varfetch = cellfetch(tetinfo{day}{epoch}(tets), 'numcells');

%returns the first tet with highest number of cells
[~,ind] = max(cell2mat(varfetch.values));

tet = tets(ind);

