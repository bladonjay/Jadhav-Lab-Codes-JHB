%This function adds cell metadata to cellinfo and spikes files for cells
%that were active on only certain epochs. Eliminates empty arrays on epochs
%where cell was not active, makes looping through cells more
%straightforward. 

%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
clear
animals = {'CS39'};
topDir = cs_setPaths();

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
    cellfilter = '$meanrate > 0';
    allcells = evaluatefilter(cellinfo,cellfilter);
    
    days = unique(allcells(:,1));
    for d = days'
        spikes = loaddatastruct(animDir, animal, 'spikes',d);
        daycells = allcells(allcells(:,1) == d,2:4);
        epochs = unique(daycells(:,1));
        daystr = getTwoDigitNumber(d);
        
        for e = epochs'
            epochcells = daycells(daycells(:,1)==e,2:3);
            
            
            for c = 1:size(epochcells)
                cell = epochcells(c,:);
                
                meanrate = spikes{d}{e}{cell(1)}{cell(2)}.meanrate;
                try
                info = spikes{d}{e}{cell(1)}{cell(2)}.info;
                catch
                    info = spikes{d}{e}{cell(1)}{cell(2)}.descript;
                end
                fields = spikes{d}{e}{cell(1)}{cell(2)}.fields;
                peak_amplitude = spikes{d}{e}{cell(1)}{cell(2)}.peak_amplitude;
                
                numspikes = cellinfo{d}{e}{cell(1)}{cell(2)}.numspikes;
                csi = cellinfo{d}{e}{cell(1)}{cell(2)}.csi;
                propbursts = cellinfo{d}{e}{cell(1)}{cell(2)}.propbursts;

                %add these to spikes and cellinfo files - just to get around
                %problem of empty cells. make data field empty in spikes, and do
                %not replace tags in cellinfo.
                
                for ep = epochs'
                    if isempty(spikes{d}{ep}{cell(1)}{cell(2)})
                        spikes{d}{ep}{cell(1)}{cell(2)}.data = [];
                        spikes{d}{ep}{cell(1)}{cell(2)}.meanrate = meanrate;
                        spikes{d}{ep}{cell(1)}{cell(2)}.info = info;
                        spikes{d}{ep}{cell(1)}{cell(2)}.fields = fields;
                        spikes{d}{ep}{cell(1)}{cell(2)}.peak_amplitude = peak_amplitude;
                    end
                    
                    if length(cellinfo{d}{ep}) < cell(1) || isempty(cellinfo{d}{ep}{cell(1)}) || ... 
                            (length(cellinfo{d}{ep}{cell(1)}) < cell(2)) || isempty(cellinfo{d}{ep}{cell(1)}{cell(2)})
                        cellinfo{d}{ep}{cell(1)}{cell(2)}.numspikes = numspikes;
                        cellinfo{d}{ep}{cell(1)}{cell(2)}.csi = csi;
                        cellinfo{d}{ep}{cell(1)}{cell(2)}.propbursts = propbursts;
                    end
                end
            end
        end
        save([animDir, animal,'spikes',daystr],'spikes');
    end
    
    save([animDir, animal,'cellinfo'],'cellinfo');
    
end
