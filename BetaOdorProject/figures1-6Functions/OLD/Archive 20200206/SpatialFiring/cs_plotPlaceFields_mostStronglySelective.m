%plot place fields for most selective cells

clear
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1'};
dataDir = [topDir,'AnalysesAcrossAnimals\'];

numcells = 5;
threshocc = 0.04;

for r = 1:length(regions)
    region = regions{r};
    
    %get list of SIs
    load([dataDir, 'selectiveCells_',region])
    %only use CS31, CS33, CS34, CS35, CS44
    selectivecells = selectivecells(ismember(selectivecells(:,1),[1 2 3 4 8]),:);
    allSI = [];
    for c = 1:size(selectivecells,1)
        cell = selectivecells(c,:);
        animal = animals{cell(1)};
        day = cell(2);
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
        try si = cellinfo{day}{2}{cell(3)}{cell(4)}.SI;
        catch ME
        end
        allSI = [allSI;si];
    end
    
    %get most left selective and most right selective
    [sorted,inds] = sort(allSI);
    cellinds = [inds(1:numcells); inds(end-(numcells-1):end)];
    sis = [sorted(1:numcells); sorted(end-(numcells-1):end)];
    
    for c = 1:length(cellinds)
        ind = cellinds(c);
        cell = selectivecells(ind,:);
        animal = animals{cell(1)};
        d = cell(2);
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        mapfields = loaddatastruct(animDir, animal, 'mapfields',d);
        
        epochs = find(~cellfun(@isempty, mapfields{d}));
        clear fields; clear occ;
        spikecount = 0;
        
        si = sis(c);
        
        for ep = 1:length(epochs)
            epoch = epochs(ep);
            if length(mapfields{d}{epoch}) >= cell(3) && length(mapfields{d}{epoch}{cell(3)}) >= cell(4) ...
                    && ~isempty(mapfields{d}{epoch}{cell(3)}{cell(4)})
                sp = mapfields{d}{epoch}{cell(3)}{cell(4)}.smoothedspikerate;
                smo =  mapfields{d}{epoch}{cell(3)}{cell(4)}.smoothedoccupancy;
                cnt = sum(mapfields{d}{epoch}{cell(3)}{cell(4)}.spikes(:));
                spikecount = spikecount + cnt;
                fields{ep} = sp;
                occ{ep} = smo;
            else
                continue
            end
        end
        occ = occ(~cellfun(@isempty,occ));
        fields = fields(~cellfun(@isempty,fields));
        [nrows, ncols] = cellfun(@size, fields);
        height = max(nrows);
        width = max(ncols);
        
        for s = 1:length(fields)
            diffheight = height-nrows(s);
            diffwidth = width-ncols(s);
            
            %pad matrices if they are different sizes before taking
            %mean. pad on all sides so that image is centered
            %if difference in size is an odd number, add one less
            %row/column at the end
            pad = [diffheight/2, diffwidth/2];
            
            %add rows
            if rem(diffheight(1), 2) == 0
                padmat = repmat(-1, pad(1), size(fields{s},2));
                if ~isempty(padmat)
                    fields{s} = [padmat; fields{s}; padmat];
                    occ{s} = [padmat; occ{s}; padmat];
                    
                end
            else
                pad(1) = round(pad(1));
                padmat1 = repmat(-1, pad(1), size(fields{s},2));
                padmat2 = repmat(-1, pad(1)-1, size(fields{s},2));
                fields{s} = [padmat1; fields{s}; padmat2];
                occ{s} = [padmat1; occ{s}; padmat2];
            end
            
            %columns
            if rem(diffwidth, 2) == 0
                padmat = repmat(-1, size(fields{s},1), pad(2));
                if ~isempty(padmat)
                    fields{s} = [padmat, fields{s}, padmat];
                    occ{s} = [padmat, occ{s}, padmat];
                end
            else
                pad(2) = round(pad(2));
                padmat1 = repmat(-1, size(fields{s},1), pad(2));
                padmat2 = repmat(-1, size(fields{s},1), pad(2)-1);
                fields{s} = [padmat1, fields{s}, padmat2];
                occ{s} = [padmat1, occ{s}, padmat2];
            end
            
        end
        
        field = mean(cat(3,fields{:}),3);
        occupancy = mean(cat(3,occ{:}),3);
        
        zero = find(occupancy <= threshocc);
        field(zero) = -1;
        
        
        imagesc(field)
        colorbar
        
        SIstr = [' SI = ',num2str(round(si,2))];
        
        titlestr = [SIstr, ' ',num2str(spikecount), ' spikes'];
        title(titlestr);
        
        figfile = [figDir, 'PlaceFields\MostSelective\',animal, '_', num2str(d), '-',num2str(cell(3)),'-',num2str(cell(4))];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        
    end
    
end