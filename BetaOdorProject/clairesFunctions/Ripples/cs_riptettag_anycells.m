function cs_riptettag_anycells(animals)
%assigns 'riptet' description to tetinfo for all tetrodes that had cells
%for each epoch

%animals = {'CS44'};
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

topDir = cs_setPaths();

for a = 1:length(animals)
    animal = animals{a};
    disp(['Tagging riptets animal ',animal]);
    dataDir = [topDir,animal,'Expt',filesep,animal,'_direct',filesep];

    load([dataDir,animal,'tetinfo.mat']);
    load([dataDir,animal,'cellinfo.mat']);
    tetfilter = '(isequal($area,''CA1''))';

    %remove old descrip fields
        filt = '~isempty($descrip)';
        old = evaluatefilter(cellinfo,filt);
        
        for f = 1:size(old,1)
            cell = old(f,:);
            if isfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)},'descrip')
                cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)} = rmfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)}, 'descrip');
            end
            if isfield(tetinfo{cell(1)}{cell(2)}{cell(3)},'descrip')
                tetinfo{cell(1)}{cell(2)}{cell(3)} = rmfield(tetinfo{cell(1)}{cell(2)}{cell(3)}, 'descrip');
            end
        end
    
    days = find(~cellfun(@isempty, cellinfo));

for d = 1:length(days)
    day = days(d);
    epochs = find(~cellfun(@isempty,cellinfo{day}));

    for ep = 1:length(epochs)
        epoch = epochs(ep);
        
        %all CA1 tets that have cells are riptets
        cellfilter = 'isequal($area,''CA1'') & (isequal($tag,''accepted'') | isequal($tag,''mua''))';
        cells = evaluatefilter(cellinfo{day}{epoch},cellfilter);
        riptets = unique(cells(:,1));
        
        for t = 1:length(riptets)
            tet = riptets(t);
            tetinfo{day}{epoch}{tet}.descrip = 'riptet';
        end
    end
end
save([dataDir,animal,'tetinfo.mat'],'tetinfo');
save([dataDir, animal,'cellinfo.mat'],'cellinfo');
end