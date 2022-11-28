clear all

%animals = {'CS31','CS33','CS34'};
animals = {'CS35'}
dataDir = 'D:\OdorPlaceAssociation\';

for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [dataDir, animal,'Expt\',animal,'_direct\'];
    
    load([animDir,animal,'cellinfo.mat'])
    load([animDir,animal,'tetinfo.mat'])
    
    masterTable = [];
    allDaysCA1 = []; allDaysPFC = [];
    
    CA1tets = evaluatefilter(tetinfo{1,1}{1,1}, 'isequal($area,''CA1'')');
    PFCtets = evaluatefilter(tetinfo{1,1}{1,1}, 'isequal($area,''PFC'')');
    
    spacer = {[]};
    
    for d = 1:length(cellinfo) %days
        
        epochs = length(cellinfo{1,d});
        
        %table = zeros(length(CA1tets)+length(PFCtets),epochs);
        load([animal,'odorTriggers0',num2str(d),'.mat'])
        odorTriggers = odorTriggers{1,d};
        
        runepochs = find(~cellfun(@isempty,odorTriggers));
        
        
        for t = 1:length(CA1tets)
            
            for ep = 1:length(runepochs)
                epoch = runepochs(ep);
                
                if length(cellinfo{1,d}{1,epoch}) >= CA1tets(t)
                    if ~isempty(cellinfo{1,d}{1,epoch}{1,CA1tets(t)})
                        cellsEachEpoch(ep) = sum(~cellfun(@isempty, cellinfo{1,d}{1,epoch}{1,CA1tets(t)}));
                    else
                        cellsEachEpoch(ep) = 0;
                    end
                else
                    cellsEachEpoch(ep) = 0;
                end
            end
            
            numcells = max(cellsEachEpoch);
            clear cellsEachEpoch
            CA1cells(t) = numcells;
            
        end
        
        for t = 1:length(PFCtets)
            
            for ep = 1:length(runepochs)
                epoch = runepochs(ep);
                    if length(cellinfo{1,d}{1,epoch}) >= PFCtets(t)
                        if ~isempty(cellinfo{1,d}{1,epoch}{1,PFCtets(t)})
                        cellsEachEpoch(ep) = sum(~cellfun(@isempty,cellinfo{1,d}{1,epoch}{1,PFCtets(t)}));
                        else
                        cellsEachEpoch(ep) = 0;
                        end
                    else 
                        cellsEachEpoch(ep) = 0;
                    end
            end
            
            numcells = max(cellsEachEpoch);
            clear cellsEachEpoch
            PFCcells(t) = numcells;
            
        end
        
        
        trials = 0;
        for r = 1:length(runepochs)
                trials = trials + length(odorTriggers{1,runepochs(r)}.allTriggers);
        end
        
        
            
        headings{1,d} = ['Day',num2str(d)];
        headings{2,d} = [num2str(length(runepochs)),' run sessions'];
        headings{3,d} = [num2str(trials),' trials'];
        
        
        
        allDaysCA1 = [allDaysCA1, CA1cells'];
        allDaysPFC = [allDaysPFC, PFCcells'];
        
    end
    
    CA1tets = CA1tets(any(allDaysCA1,2));
    allDaysCA1 = allDaysCA1(any(allDaysCA1,2),:);
    
    PFCtets = PFCtets(any(allDaysPFC,2));
    allDaysPFC = allDaysPFC(any(allDaysPFC,2),:);
        
    totalCA1 = sum(sum(allDaysCA1));
    totalPFC = sum(sum(allDaysPFC));
    
    masterTable = [num2cell(allDaysCA1); repmat(spacer,1,d); num2cell(allDaysPFC)]; 
    
    
    tetlabels = [spacer;spacer;{'CA1tetrodes'};num2cell(CA1tets);{'PFCtetrodes'};num2cell(PFCtets)];
    
    masterTable = [headings; masterTable];
    masterTable = [tetlabels, masterTable];
    
    masterTable{(length(CA1tets)+3), d+2} = 'Total CA1 cells = ';
    masterTable{(length(CA1tets)+3), d+3} = totalCA1;
    
    masterTable{(length(CA1tets)+4 + length(PFCtets)), d+2} = 'Total PFC cells = ';
    masterTable{(length(CA1tets)+4 + length(PFCtets)), d+3} = totalPFC;
    
    
    xlswrite([dataDir,'CellMetaData.xlsx'], masterTable ,animal);
    
    clear CA1cells; clear allDaysCA1; clear PFCcells; clear allDaysPFC;
end