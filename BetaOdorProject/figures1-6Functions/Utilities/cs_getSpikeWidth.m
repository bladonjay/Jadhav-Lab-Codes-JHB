clear
topDir = cs_setPaths;
animals = {'CS44'};
method = 2;
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
    cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
    waves = loaddatastruct(animDir, animal,'waves');
    spikes = loaddatastruct(animDir, animal,'spikes');
    
    days = find(~cellfun(@(x) isempty(x), waves));
    
    
    for day = days
        
        daystr = getTwoDigitNumber(day);
        %waves = loaddatastruct(animDir, animal, 'waves',day);
        disp(['Doing ',animal,' day ',daystr]);
        
        eps = find(~cellfun(@(x) isempty(x), spikes{day}));
        for ep = eps
            
            tets = find(~cellfun(@(x) isempty(x), spikes{day}{ep}));
            
            for tet = tets
                
                cells = find(~cellfun(@(x) isempty(x), spikes{day}{ep}{tet}));
                
                for cell = cells
                    curWaves = waves{day}{1}{tet}{cell};
                    
                    [peak,peakIdx]=(max(curWaves,[],2)); %get peak for each spike
                    [trou,trouIdx]=(min(curWaves,[],2)); %get trough for each spike
                    
                    if method ==1 
                    %method 1: width between peak and trough
                    diff = abs(trouIdx-peakIdx);
                    spikewidth= mean(diff*(1/30000)*1000); % in ms
                    end
                    
                    if method ==2
                    %method 2: width at half height
                    meanwave = mean(curWaves);
                    upsamp   = 10;
                    dummy    = interp1(meanwave,linspace(1,length(meanwave),length(meanwave)*upsamp));
                     [peak, indpeak] = min(dummy);
                     first    = find(dummy<peak/2,1,'first');
                     last     = find(dummy(first+1:end)>peak/2,1,'first')+first;
                     diff      = last-first;
                     spikewidth= mean(diff*(1/30000)*1000); % in ms
%                      
%                      plot(dummy),hold on
%                      plot([first first],[-10^3 10^3],'g:')
%                      plot([last last],[-10^3 10^3],'g:')
%                      plot([0 numel(dummy)],[peak/2 peak/2],'g:')
%                      plot([first last],[peak/2 peak/2],'g')
                    end
                      
                    
                    
                    
                    spikes{day}{ep}{tet}{cell}.spikewidth = spikewidth;
                    cellinfo{day}{ep}{tet}{cell}.spikewidth = spikewidth;
                end
            end
        end
        
        save([animDir, animal,'spikes',daystr],'spikes');
    end
    save([animDir, animal,'cellinfo'],'cellinfo');
end
