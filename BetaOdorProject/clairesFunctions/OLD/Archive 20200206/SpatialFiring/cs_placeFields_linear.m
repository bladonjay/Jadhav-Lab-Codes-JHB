%cs_placeFields_linear

%Get linear place fields for right and left trajectories. Do for selective
%and non-selective cells. Subtract the two to get the spatial selectivity
%index for each spatial bin. Indicate on plots the NP, turn away from NP,
%and choice point. 

%traj 1 = left traj, traj 3 = right traj. 

%ADD: look at correct trials only. anneal run time with NP time (speed
%filter is messing this up).
%Then compare to incorrect trials. 

%Also look only at reward zone


clear

speedthresh = 5;
binsize = 1; % cm square 
%stdv = 3; 
%threshocc = 0.02; % Threshold occupancy in seconds
%npExcludeDist = 30; %cm away from NP


animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC'};
topDir = cs_setPaths;

riptetfilter = '(isequal($descrip, ''riptet''))';
%timefilter = { {'DFTFsj_getlinstate', '(($state ~= 0) & (abs($linearvel) >= 5))', 3, 'headdir',1}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',1} };


for r = 1:length(regions)
    region = regions{r};
    load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region,'.mat']);
    %load([topDir,'AnalysesAcrossAnimals\npCells_',region,'.mat']);
    
    SIs = [];
    trajcomps = [];
    for a = 1:length(animals)
        animal = animals{a};
        animcells = selectivecells(selectivecells(:,1) == a, 2:4);
        %animcells = npCells(npCells(:,1) == a, 2:4);
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        load([animDir,animal,'cellinfo.mat']);
        days = unique(animcells(:,1));
        
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            load([animDir,animal,'spikes',daystr,'.mat']);
            load([animDir,animal,'pos',daystr,'.mat']);
            load([animDir,animal,'ripples',daystr,'.mat']);
            load([animDir,animal,'linpos',daystr,'.mat']);
            epochs = find(~cellfun(@isempty,linpos{day}));
            runmatrix = [repmat(day, length(epochs),1), epochs'];
            
            linstate = DFTFsj_getlinstate(animDir,animal,runmatrix, 6);        
            rips = DFTFsj_getriptimes(animDir,animal,runmatrix, 'tetfilter',riptetfilter,'minthresh',3);
            
            
            daycells = animcells(animcells(:,1) == day,[2,3]);
            
            clear newpos
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                
                excludePeriods = [];
                tmp = pos{day}{epoch}.data(:,1); %time series, exclude times from here
                tmp_linpos = linpos{day}{epoch}.statematrix.linearDistanceToWells(:,1); %linear distance from NP
                postime = linpos{day}{epoch}.statematrix.time;
                %posdata = pos{day}{epoch}.data(:,[2,3]);
                timestep = tmp(2) - tmp(1);
                % excludePeriods = getExcludePeriods(time, included)
                %  Calculates the start and end times for all exclusion periods, given a
                % time vector and an include vector of 1s and 0s of the same length.
                
                
                excludeSpeed = getExcludePeriods(postime,abs(linstate{day}{epoch}.linearvel) > speedthresh);
                tmp = tmp(~isExcluded(tmp, excludeSpeed));
%                 
%                
%                 npCoord = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
%                 
%                 excludePos = getExcludePeriods(postime,posdata(:,2) < (npCoord(1,2) - npExcludeDist));
%                 tmp = tmp(~isExcluded(tmp, excludePos));

                excludeRipples = getExcludePeriods(rips{day}{epoch}.time', double(~rips{day}{epoch}.nripples));
                
                tmp = tmp(~isExcluded(tmp, excludeRipples));
                
                %left trajectories
                excludeState_Left = getExcludePeriods(postime,(linstate{day}{epoch}.state == 1 ));
                tmp_l = tmp(~isExcluded(tmp, excludeState_Left));
                
                %right trajectories
                excludeState_Right = getExcludePeriods(postime,(linstate{day}{epoch}.state == 3 ));
                tmp_r = tmp(~isExcluded(tmp, excludeState_Right));
 
                posind_l = lookup(tmp_l, postime);
                posind_r = lookup(tmp_r, postime);
                
%                 goodpos_l = posdata(posind_l,:);
%                 goodpos_r = posdata(posind_r,:);
                goodpos_l = tmp_linpos(posind_l);
                goodpos_r = tmp_linpos(posind_r);
                
                %figure, plot(goodpos(:,1), goodpos(:,2),'k.');
                
%                 state = find(linstate{day}{epoch}.state > 0);
%                 posinclude = intersect(ls, state);
%                 goodpos = pos{day}{epoch}.data(posinclude,:);  
                
                newpos_l{epoch} = goodpos_l;
                newpos_r{epoch} = goodpos_r;
                
                excl{epoch}.rip = excludeRipples;
                excl{epoch}.sp = excludeSpeed;
                %excl{epoch}.pos = excludePos;
                excl{epoch}.left = excludeState_Left;
                excl{epoch}.right = excludeState_Right;
                
%                 figure, plot(goodpos_l(:,1), goodpos_l(:,2),'k.');
%                 figure, plot(goodpos_r(:,1), goodpos_r(:,2),'k.');
%                 figure, plot(pos{day}{epoch}.data(:,2), pos{day}{epoch}.data(:,3),'k.');
            end

            for c = 1:size(daycells,1)
                flag = 0;
                cell = daycells(c,:);
                allpos_l = []; allpos_r = [];
                allspikes_l = []; allspikes_r = [];
                si = [];
                for ep = 1:length(epochs)
                    epoch = epochs(ep);
                    if isa(spikes{day}{epoch}{cell(1)},'cell')
                        if isa(spikes{day}{epoch}{cell(1)}{cell(2)},'struct')
                            
                            if isfield(cellinfo{day}{epoch}{cell(1)}{cell(2)},'SI')
                            si = cellinfo{day}{epoch}{cell(1)}{cell(2)}.SI;
                            else
                                flag = 1;
                                continue
                            end
                            
                            tmp = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                            
                            %for each epoch, the times to exclude were
                            %already found. So use these to exclude spikes
                            %that fall within Exclude periods
                            excludeSpeed = excl{epoch}.sp;
                            tmp = tmp(~isExcluded(tmp, excludeSpeed));
                            
                            excludeRipples = excl{epoch}.rip;
                            tmp = tmp(~isExcluded(tmp, excludeRipples));
                            
%                             excludePos = excl{epoch}.pos;
%                             tmp = tmp(~isExcluded(tmp, excludePos));
                            
                            excludeState_Left = excl{epoch}.left;
                            tmp_l = tmp(~isExcluded(tmp, excludeState_Left));
                            
                            excludeState_Right = excl{epoch}.right;
                            tmp_r = tmp(~isExcluded(tmp, excludeState_Right));
                            
                            postime = linpos{day}{epoch}.statematrix.time;
                            
                            spikeind_l = lookup(tmp_l, postime);
                            %numspikes_l = length(spikeind_l);
                            spikepos_l = linpos{day}{epoch}.statematrix.linearDistanceToWells(spikeind_l,1);
                            
                            spikeind_r = lookup(tmp_r, postime);
                            %numspikes_r = length(spikeind_r);
                            spikepos_r = linpos{day}{epoch}.statematrix.linearDistanceToWells(spikeind_r,1);
                            
                            allpos_l = [allpos_l; newpos_l{epoch}];
                            allpos_r = [allpos_r; newpos_r{epoch}];
                            allspikes_l = [allspikes_l; spikepos_l];
                            allspikes_r = [allspikes_r; spikepos_r];
                            %allspikes = [allspikes; numspikes_1, numspikes_r];
                           
                            
                            
                        end
                    end
                end
                if flag == 1
                    continue
                end
                %meanspikerate = mean(allspikes,1) ./ timestep;
                
                %LEFT
                Min = floor(min(allpos_l)); 
                Max = ceil(max(allpos_l));
                bins = (Min:binsize:Max);
               
                
                [occupancy, ticks] = histcounts(allpos_l, bins);
                
                 if ~isempty(allspikes_l)
                    %1) Calculate occupancy-normalized spikerate
                    [s, BX] = histcounts(allspikes_l,  bins);
                    
                    %make sure sizes of occupancy and spike matricies are
                    %the same
%                     if length(s) ~= length(occupancy)
%                         numrows = min([length(occupancy), length(s)]);
%                         s = s(1:numrows,:);
%                         occupancy = occupancy(1:numrows,:);
%                     end
%                     
%                     if size(s,2) ~= size(occupancy,2)
%                         numcols = min([size(occupancy,2), size(s,2)]);
%                         s = s(:,1:numcols);
%                         occupancy = occupancy(:,1:numcols);
%                     end
                    
                    nonzero = find(occupancy ~= 0);
                    spikerate = zeros(size(s));
                    spikerate(nonzero) = s(nonzero) ./(occupancy(nonzero) * timestep );
                    
                    spikerate_l = smoothdata(spikerate);%spikerate = spikerate/timestep;
                     %bad = (spikerate > 100);
                    %spikerate(bad) = 0;
                    
%                     occupancy = timestep*occupancy;
%                     z = find(occupancy <= threshocc);
%                     spikerate = spikerate(~z);
%                     spikerate_l = mean(spikerate);
                    
                 else
                     spikerate_l = zeros(1,length(bins-1));
                 end
                 
                 
                 %RIGHT
                Min = floor(min(allpos_l)); 
                Max = ceil(max(allpos_l));
                bins = (Min:binsize:Max);
                
                 [occupancy, ticks] = histcounts(allpos_r, bins);
                
                if ~isempty(allspikes_r)
                    %1) Calculate occupancy-normalized spikerate
                    [s, BX] = histcounts(allspikes_r, bins);
                    
                    %make sure sizes of occupancy and spike matricies are
                    %the same
%                     if size(s,1) ~= size(occupancy,1)
%                         numrows = min([size(occupancy,1), size(s,1)]);
%                         s = s(1:numrows,:);
%                         occupancy = occupancy(1:numrows,:);
%                     end
                    
%                     if size(s,2) ~= size(occupancy,2)
%                         numcols = min([size(occupancy,2), size(s,2)]);
%                         s = s(:,1:numcols);
%                         occupancy = occupancy(:,1:numcols);
%                     end
                    
                    
                    nonzero = find(occupancy ~= 0);
                    spikerate = zeros(size(s));
                    spikerate(nonzero) = s(nonzero) ./(timestep* occupancy(nonzero) );
                    spikerate_r = smoothdata(spikerate);
%                      bad = (spikerate > 100);
%                     spikerate(bad) = 0;
%                     
%                     occupancy = timestep*occupancy;
%                     z = find(occupancy <= threshocc);
%                     spikerate = spikerate(~spikerate(z));
%                     spikerate_r = mean(spikerate);
                else 
                    spikerate_r = zeros(1,length(bins-1));;
                end
                
%                     trajcomp = (spikerate_l - spikerate_r)/(spikerate_l + spikerate_r);
%                     
%                     SIs = [SIs; si];
%                     trajcomps = [trajcomps; trajcomp];
                %indexval = (spikerate_l - spikerate_r)./(spikerate_l + spikerate_r);
                figure
                f = gcf;
                f.Position = [2000 300 800 250];
                
                
                ax1 = subplot(2,1,1);
                imagesc(spikerate_l);
                colorbar
                %cmap = [rgb('DodgerBlue'); rgb('Black'); rgb('DeepPink')];
                %ax.Colormap = [rgb('DodgerBlue'); rgb('Black'); rgb('DeepPink')];
                [cmap]=buildcmap('kc'); 
                colormap(ax1,cmap)
                
                ax2 = subplot(2,1,2);
                
                imagesc(spikerate_r);
                [cmap]=buildcmap('km'); 
                colormap(ax2, cmap)
                colorbar
                
                suptitle(['SI = ',num2str(si)]);
                pause
                close all

            end
            
        end
        
    end
    

end