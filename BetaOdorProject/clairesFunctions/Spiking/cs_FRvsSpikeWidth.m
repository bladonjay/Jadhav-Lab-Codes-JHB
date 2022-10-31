%plot FR vs Spikewidth
clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

topDir = cs_setPaths;
%days = [1:5];
regions = {'CA1','PFC'};

for r = 1:length(regions)
    region = regions{r};
    
    pyrdat = [];
    intdat = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
        
        runeps = cs_getRunEpochs(animDir, animal,'odorplace');
        days = unique(runeps(:,1));
        
        %pyr
        cellfilt = ['strcmp($area,''',region,''') &strcmp($type,''pyr'') & $meanrate > 0 & $spikewidth > 0'];
        cells = evaluatefilter(cellinfo,cellfilt);
        cells = cells(ismember(cells(:,[1 2]),runeps,'rows'),:);
        cells = unique(cells(:,[1 3 4]),'rows');
        
        
        for c = 1:size(cells)
            cell = cells(c,:);
            eps = cs_findGoodEpochs(cellinfo{cell(1)},{'meanrate','spikewidth'},cell(2:3));
            
            fr = [];
            for ep = eps'
                f = cellinfo{cell(1)}{ep}{cell(2)}{cell(3)}.meanrate;
                fr = [fr;f];
            end
            %fr = cellinfo{cell(1)}{eps(1)}{cell(2)}{cell(3)}.meanrate;
            fr = mean(fr);
            sw = cellinfo{cell(1)}{eps(1)}{cell(2)}{cell(3)}.spikewidth;
            dat = [fr, sw];
            pyrdat = [pyrdat;dat];
        end
        
        %int
        cellfilt = ['strcmp($area,''',region,''') &strcmp($type,''int'') & $meanrate > 0 & $spikewidth > 0'];
        cells = evaluatefilter(cellinfo,cellfilt);
        cells = cells(ismember(cells(:,[1 2]),runeps,'rows'),:);
        cells = unique(cells(:,[1 3 4]),'rows');
        
        
        for c = 1:size(cells)
            cell = cells(c,:);
            eps = cs_findGoodEpochs(cellinfo{cell(1)},{'meanrate','spikewidth'},cell(2:3));
            
            fr = [];
            for ep = eps'
                f = cellinfo{cell(1)}{ep}{cell(2)}{cell(3)}.meanrate;
                fr = [fr;f];
            end
            %fr = cellinfo{cell(1)}{eps(1)}{cell(2)}{cell(3)}.meanrate;
            fr = mean(fr);
            sw = cellinfo{cell(1)}{eps(1)}{cell(2)}{cell(3)}.spikewidth;
            dat = [fr, sw];
            if fr < 13
                disp(check)
            end
            intdat = [intdat;dat];
        end
    end
    figure,
        plot(pyrdat(:,1),pyrdat(:,2),'r.');
        hold on
        plot(intdat(:,1),intdat(:,2),'b.');
        set(gca, 'XScale', 'log')
        title(region)
end