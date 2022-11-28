function celleps = cs_getCellEpochs(cellinfostruct,day, tet,cell,runeps)

celleps = [];
for r = 1:length(runeps)
    runep = runeps(r);
    if length(cellinfostruct{day}{runep}{tet}) > cell
        test = cellinfostruct{day}{runep}{tet}{cell};
        if isfield(test,'numspikes')
            celleps = [celleps;runep];
        end
    end
end
