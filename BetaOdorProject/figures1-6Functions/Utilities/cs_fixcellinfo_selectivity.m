%cs_fixcellinfo

clear
topDir = cs_setPaths;

animals = {'CS31','CS33','CS34','CS35'};

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    
    load([animDir,animal,'cellinfo.mat'])
    
    filt = ['(~isempty($SI))'];
    test = evaluatefilter(cellinfo,filt);
    for f = 1:size(test,1)
        cell = test(f,:);
        if isfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)},'SI')
            cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)} = rmfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)}, 'SI');
        end
    end
    
    filt = '(strcmp($selectivity, ''leftSelective'')) || (strcmp($selectivity, ''rightSelective''))';
    test = evaluatefilter(cellinfo,filt);
    for f = 1:size(test,1)
        cell = test(f,:);
        if isfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)},'selectivity')
            cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)} = rmfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)}, 'selectivity');
        end
    end
    
    filt = ['(~isempty($SI) || (strcmp($selectivity, ''leftSelective'')) || (strcmp($selectivity, ''rightSelective'')))'];
    check = evaluatefilter(cellinfo,filt);
    
    num = size(check,1)
    
    save([topDir,animal,'Expt\',animal,'_direct\',animal,'cellinfo.mat'],'cellinfo');
    
end