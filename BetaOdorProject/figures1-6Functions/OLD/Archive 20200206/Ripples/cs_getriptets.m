function riptets = cs_getriptets(animal)

topDir = cs_setPaths();
animDir = [topDir, animal,'Expt',filesep,animal,'_direct',filesep];
tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
tetfilter = '(isequal($descrip, ''riptet''))';
riptets = evaluatefilter(tetinfo, tetfilter);
riptets = unique(riptets(:,3));
