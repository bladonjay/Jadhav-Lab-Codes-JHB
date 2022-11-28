%classify place cells
%it is a place cell if it has a peak rate of >3Hz on track (maybe make it
%5Hz to be more conservative?)

%Can tag cellinfo

%List selective cells that are place cells
clear; close all
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','','','','CS44'};
threshold = 3;

for a = [1 2 3 4 8]
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    mapfiles = dir([animDir,animal,'mapfields*']);
    cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
    %Remove existing tags
    filt2 = '~isempty($placetag)';
    old = evaluatefilter(cellinfo,filt2);
    
    for f = 1:size(old,1)
        cell = old(f,:);
        if isfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)},'placetag')
            cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)} = rmfield(cellinfo{cell(1)}{cell(2)}{cell(3)}{cell(4)}, 'placetag');
        end
    end
    
    
    for d = 1:length(mapfiles)
        load([animDir,mapfiles(d).name]);
        day = length(mapfields);
        epochs = find(~cellfun(@isempty, mapfields{d}));
        
        filt = 'length($smoothedspikerate) > 0';
        cells = evaluatefilter(mapfields{d},filt);
        cells = unique(cells(:,2:3),'rows');
        
        
        
        for c = 1:size(cells,1)
            cell = cells(c,:);
            peaks = [];
            goodeps = [];
            for ep = epochs
                if length(mapfields{d}{ep})< cell(1) || length(mapfields{d}{ep}{cell(1)})< cell(2) ...
                        || isempty(mapfields{d}{ep}{cell(1)}) || isempty(mapfields{d}{ep}{cell(1)}{cell(2)})
                    continue
                else
                    
                    field = mapfields{d}{ep}{cell(1)}{cell(2)}.smoothedspikerate;
                    numspikes = sum(mapfields{d}{ep}{cell(1)}{cell(2)}.spikes(:));
                    
                    peak = max(field(:));
                    peaks = [peaks, peak];
                    goodeps = [goodeps, ep];
                end
                
            end
            
            chk = find(peaks >= threshold);
            
            for j = goodeps
                if any(chk)
                    cellinfo{d}{j}{cell(1)}{cell(2)}.placetag = 'placecell';
                else
                    cellinfo{d}{j}{cell(1)}{cell(2)}.placetag = 'nonresp';
                end
            end
            
        end
        daystr = getTwoDigitNumber(day);
        save([animDir, animal, 'cellinfo'],'cellinfo');
        
    end
end

cs_listPlaceCells;