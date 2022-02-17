%cs_placeFields

%next, compare avg fr between right and left traj, get a value for this.
%then compare to odor selectivity index. Are they correlated? 
%Can include only selective cells, or both selective and non-selective cells...
%traj 1 = left traj, traj 3 = right traj. limit to these. 
%add a position filter to exclude area around NP, often cells will prefer
%one turn over another, don't want to include this, just well trajectories.
%


clear

speedthresh = 20;
binsize = 1; % cm square 
std = 3; 
threshocc = 0.02; % Threshold occupancy in seconds
npExcludeDist = 30; %cm square around NP


animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'}; 
animNums = [1,2,3,4,8]; %don't use CS39, CS41, CS42, not enough pos (short track)
regions = {'CA1'};
topDir = cs_setPaths;

riptetfilter = '(isequal($descrip, ''riptet''))';
%timefilter = { {'DFTFsj_getlinstate', '(($state ~= 0) & (abs($linearvel) >= 5))', 3, 'headdir',1}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',1} };


for r = 1:length(regions)
    region = regions{r};
    load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region,'.mat']);
    
    %load([topDir,'AnalysesAcrossAnimals\npCells_',region,'.mat']);

    for a = animNums
        animal = animals{a};
        animcells = selectivecells(selectivecells(:,1) == a, 2:4);
        %animcells = npCells(npCells(:,1) == a, 2:4);
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        
        load([animDir, animal,'cellinfo.mat']);
        load([animDir,animal,'rippletimes.mat']);
        days = unique(animcells(:,1));
        
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            load([animDir,animal,'spikes',daystr,'.mat']);
            load([animDir,animal,'pos',daystr,'.mat']);
            load([animDir,animal,'linpos',daystr,'.mat']);
            epochs = find(~cellfun(@isempty,linpos{day}));
            runmatrix = [repmat(day, length(epochs),1), epochs'];
            
            linstate = DFTFsj_getlinstate(animDir,animal,runmatrix, 6);        
            %rips = DFTFsj_getriptimes(animDir,animal,runmatrix, 'tetfilter',riptetfilter,'minthresh',3);
            
            daycells = animcells(animcells(:,1) == day,[2,3]);
            
            clear newpos
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                
                rips = [ripple{day}{epoch}.starttime, ripple{day}{epoch}.endtime];
                
                excludePeriods = [];
                tmp = pos{day}{epoch}.data(:,1);
                postime = pos{day}{epoch}.data(:,1);
                posdata = pos{day}{epoch}.data(:,[2,3]);
                timestep = tmp(2) - tmp(1);
                % excludePeriods = getExcludePeriods(time, included)
                %  Calculates the start and end times for all exclusion periods, given a
                % time vector and an include vector of 1s and 0s of the same length.
                excludeSpeed = getExcludePeriods(postime,abs(linstate{day}{epoch}.linearvel) > speedthresh);
                tmp = tmp(~isExcluded(tmp, excludeSpeed));
%                 
%                 excludeState = getExcludePeriods(pos{day}{epoch}.data(:,1),linstate{day}{epoch}.state > 0);
%                 tmp = tmp(~isExcluded(tmp, excludeState));
                npCoord = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
                %npbox = [npCoord(1,1)-npExcludeDist, npCoord(1,2)- 2*npExcludeDist; npCoord(1,1)+npExcludeDist, npCoord(1,2)+ npExcludeDist];
                %1st column x values, 2nd column y values
                
                %((posdata(:,1) < npbox(1,1) | posdata(:,1) > npbox(2,1)) & ... 
                excludePos = getExcludePeriods(postime,posdata(:,2) < (npCoord(1,2) - npExcludeDist));
                tmp = tmp(~isExcluded(tmp, excludePos));


                excludeState = getExcludePeriods(postime,(linstate{day}{epoch}.state == 1 | linstate{day}{epoch}.state == 3));
                tmp = tmp(~isExcluded(tmp, excludeState));
%                 
                %excludeRipples = getExcludePeriods(rips{day}{epoch}.time', double(~rips{day}{epoch}.nripples));
                tmp = tmp(~isExcluded(tmp, rips));
                
                
 
                posind = lookup(tmp, postime);
                goodpos = pos{day}{epoch}.data(posind,[2,3]);
                %figure, plot(goodpos(:,1), goodpos(:,2),'k.');
                
%                 state = find(linstate{day}{epoch}.state > 0);
%                 posinclude = intersect(ls, state);
%                 goodpos = pos{day}{epoch}.data(posinclude,:);  
                
                newpos{epoch} = goodpos;
                
%                 figure, plot(newpos(:,2), newpos(:,3),'k.');
%                 figure, plot(pos{day}{epoch}.data(:,2), pos{day}{epoch}.data(:,3),'k.');
            end

            for c = 1:size(daycells,1)
                cell = daycells(c,:);
                allpos = [];
                allspikes = []; %plot one figure for whole day
                for ep = 1:length(epochs)
                    epoch = epochs(ep);
                    if isa(spikes{day}{epoch}{cell(1)},'cell')
                        if isa(spikes{day}{epoch}{cell(1)}{cell(2)},'struct')
                            
                            if isfield(cellinfo{day}{epoch}{cell(1)}{cell(2)},'SI')
                            SI = cellinfo{day}{epoch}{cell(1)}{cell(2)}.SI;
                            selectivity = cellinfo{day}{epoch}{cell(1)}{cell(2)}.selectivity;
                            else
                                SI = '';
                            end
                            
                            tmp = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                            
                            %for each epoch, the times to exclude were
                            %already found. So use these to exclude spikes
                            %that fall within Exclude periods
                            tmp = tmp(~isExcluded(tmp, excludeSpeed));
                            %tmp = tmp(~isExcluded(tmp, excludePos));
                            tmp = tmp(~isExcluded(tmp, excludeState));
                            tmp = tmp(~isExcluded(tmp, rips));
                            
                            spikeind = lookup(tmp, spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1));
                            spikepos = spikes{day}{epoch}{cell(1)}{cell(2)}.data(spikeind,[2,3]);
                            
                            allpos = [allpos; newpos{epoch}];
                            allspikes = [allspikes; spikepos];
                            
                        end
                    end
                end
                
                %should this be calculated separately for each cell? Should
                %be the same in most cases, except for when a cell only
                %spiked in one run epoch and not another  (i.e. drift)
                minx = floor(min(allpos(:,1))) - 1; 
                maxx = ceil(max(allpos(:,1))) + 1;
                binx = (minx:binsize:maxx);
                miny = floor(min(allpos(:,2))) - 1;
                maxy = ceil(max(allpos(:,2))) + 1;
                biny = (miny:binsize:maxy);
                
                [occupancy, xticks, yticks] = hist2(allpos(:,1), allpos(:,2), binx, biny);
                
                if ~isempty(allspikes)
                    
                    
 %                  %1) Calculate occupancy-normalized spikerate
                    [s, BX, BY] = hist2(allspikes(:,1), allspikes(:,2), binx, biny);
                    
                    %make sure sizes of occupancy and spike matricies are
                    %the same
                    if size(s,1) ~= size(occupancy,1)
                        numrows = min([size(occupancy,1), size(s,1)]);
                        s = s(1:numrows,:);
                        occupancy = occupancy(1:numrows,:);
                    end
                    
                    if size(s,2) ~= size(occupancy,2)
                        numcols = min([size(occupancy,2), size(s,2)]);
                        s = s(:,1:numcols);
                        occupancy = occupancy(:,1:numcols);
                    end
                    
                    
                    nonzero = find(occupancy ~= 0);
                    spikerate = zeros(size(s));
                    spikerate(nonzero) = s(nonzero) ./(timestep* occupancy(nonzero) );
                    
                    bad = find(spikerate > 100);
                    spikerate(bad) = 0;

 %                  % 2) Smooth spikerate and occupancy separately
                    g = gaussian2(std,(5*std)); % three pixel kernel... thats huge
                    smoothedspikerate = filter2(g,(spikerate)); % is this the right filter? yeah but its oversmoothing
                    %smoothedoccupancy = [];
                    %smoothedoccupancy = zeros(size(output.spikes));
                    smoothedoccupancy = filter2(g, occupancy);
                    
 %                  % 3) Turn occupancy to seconds and set spikerate wherever occupancy
                    %is < threshold occupancy in seconds to 0
                    
                    occupancy = timestep*occupancy;
                    smoothedoccupancy = timestep*smoothedoccupancy;
                    
                    %zero = find(smoothedoccupancy == 0);
                    % zero out places he hasnt visited enough
                    zout = find(smoothedoccupancy <= threshocc);
                    %zz = find(smoothedspikerate <= 0);
                    %output.smoothedspikerate(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
                    smoothedoccupancy(zout) = -1;
                    
                    smoothedspikerate(zout) = -1;
                    
                    
                    cmax = max(smoothedspikerate(:));
                    
                    imagesc(smoothedspikerate, [-1, cmax])
                    colorbar
                    title([animal, ' ', region, ' ', num2str(day), ' ', num2str(cell(1)), ' ', num2str(cell(2)), ' SI= ', num2str(SI), ' ', selectivity]); 
                    pause
                    close
                    

                end
                    
            end
            
        end
        
    end
end