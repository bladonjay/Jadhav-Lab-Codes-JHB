%cs_plotPlaceFields
clear
[topDir, figDir] = cs_setPaths();
dataDir = [topDir, 'AnalysesAcrossAnimals\'];

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
region = 'CA1';
threshocc = 0.04; % Threshold occupancy in seconds

%plot all place cells, or plot only selective cells
%cellstoplot = 'selectiveCells';
cellstoplot = 'placecells';

switch cellstoplot
    case 'selectiveCells'
        load([dataDir,'selectiveCells_',region]);
        cells = selectivecells;
    case 'placecells'
        load([dataDir,'placecells_',region]);
        cells = placecells;
end

for a = [1 2 3 4 8]
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
    animCells = cells(cells(:,1) == a, 2:4);
    days = unique(animCells(:,1));
    for d = days'
        mapfields = loaddatastruct(animDir, animal, 'mapfields',d);
        daycells = animCells(animCells(:,1)==d,2:3);
        for c = 1:size(daycells)
            cell = daycells(c,:);
            epochs = find(~cellfun(@isempty, mapfields{d}));
            clear fields; clear occ;
            spikecount = 0;
            
            if strcmp(cellstoplot, 'selectiveCells')
                si = cellinfo{d}{epochs(1)}{cell(1)}{cell(2)}.SI;
            else
                si = [];
            end
            fields = {};
            occ = {};
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                if length(mapfields{d}{epoch}) >= cell(1) && length(mapfields{d}{epoch}{cell(1)}) >= cell(2) ...
                        && ~isempty(mapfields{d}{epoch}{cell(1)}{cell(2)})
                    sp = mapfields{d}{epoch}{cell(1)}{cell(2)}.smoothedspikerate;
                    smo =  mapfields{d}{epoch}{cell(1)}{cell(2)}.smoothedoccupancy;
                    cnt = sum(mapfields{d}{epoch}{cell(1)}{cell(2)}.spikes(:));
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
            if ~isempty(si)
                SIstr = [' SI = ',num2str(round(si,2))];
            else
                SIstr ='';
            end
            

            cellnum = find(ismember(cells,[a,d,cell(1),cell(2)],'rows'));
            
            if strcmp(cellstoplot, 'placecells')
            placemaps.maps{cellnum,1} = field;
            placemaps.inds{cellnum,1} = [a,d,cell(1),cell(2)];
            end
            
            titlestr = ['Cell ', num2str(cellnum), SIstr, ' ',num2str(spikecount), ' spikes'];
            title(titlestr);
            
            figfile = [figDir, 'PlaceFields\',cellstoplot,'\',animal, '_', num2str(d), '-',num2str(cell(1)),'-',num2str(cell(2))];
            
            savefig(figfile)
            print('-djpeg', figfile);
            print('-dpdf', figfile);
        end
    end
end

if exist('placemaps')
    save([dataDir, 'placemaps'],'placemaps')
end

