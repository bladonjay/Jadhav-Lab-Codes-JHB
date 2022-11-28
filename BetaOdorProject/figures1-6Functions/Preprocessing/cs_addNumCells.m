function cs_addNumCells(animdirect,fileprefix)

try 
    load([animdirect, fileprefix,'cellinfo']);
    cellfilter = '($numspikes > 0)';
    cells = evaluatefilter(cellinfo,cellfilter);
catch
    warning('cellinfo structure does not exist')
end

try
    load([animdirect, fileprefix, 'tetinfo']);
    tetfilter = '(~isempty($numcells))';
    tets = evaluatefilter(tetinfo,tetfilter);
    for t = 1:size(tets,1)
        tetcells = sum(ismember(cells(:,1:3),tets(t,:),'rows'));
        tetinfo{tets(t,1)}{tets(t,2)}{tets(t,3)}.numcells = tetcells;
    end
     save([animdirect, fileprefix,'tetinfo'], 'tetinfo');
catch
    warning('tetinfo structure does not exist')
end

