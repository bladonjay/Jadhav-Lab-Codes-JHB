%cs_placeFields_trajcompare

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
npExcludeDist = 30; %cm away from NP


animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
topDir = cs_setPaths;

riptetfilter = '(isequal($descrip, ''riptet''))';
timefilter = { {'DFTFsj_getlinstate', '(($state ~= 0) & (abs($linearvel) >= 5))', 3, 'headdir',1}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',1} };


for r = 1:length(regions)
    region = regions{r};
    %load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region,'.mat']);
    load([topDir,'AnalysesAcrossAnimals\npCells_',region,'.mat']);
    
    SIs = [];
    trajcomps = [];
    for a = [1 2 3 4 8] %1:length(animals)  Don't use animals that had truncated track
        animal = animals{a};
        %animcells = selectivecells(selectivecells(:,1) == a, 2:4);
        animcells = npCells(npCells(:,1) == a, 2:4);
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
            runmatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
            load([animDir,animal,'DIO',daystr,'.mat']);
            
            runmatrix = runmatrix(runmatrix(:,1) == day,:);
            epochs = runmatrix(:,2);

            
            linstate = DFTFsj_getlinstate(animDir,animal,runmatrix, 6);        
            rips = DFTFsj_getriptimes(animDir,animal,runmatrix, 'tetfilter',riptetfilter,'minthresh',3);
            
            
            daycells = animcells(animcells(:,1) == day,[2,3]);
            
            clear newpos
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                
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
%                
                npCoord = linpos{day}{epoch}.wellSegmentInfo.wellCoord;
                
                excludePos = getExcludePeriods(postime,posdata(:,2) < (npCoord(1,2) - npExcludeDist));
                tmp = tmp(~isExcluded(tmp, excludePos));
                
                excludeReward_l = cs_parseDIOtimes(animal, day, epoch, 1);
                excludeReward_r = cs_parseDIOtimes(animal, day, epoch, 2);
                excludeReward = sort([excludeReward_l; excludeReward_r],1);
                tmp = tmp(~isExcluded(tmp, excludeReward));
%                
                excludeRipples = getExcludePeriods(rips{day}{epoch}.time', double(~rips{day}{epoch}.nripples));
                tmp = tmp(~isExcluded(tmp, excludeRipples));
                
                excludeState_Left = getExcludePeriods(postime,(linstate{day}{epoch}.state == 1 ));
                tmp_l = tmp(~isExcluded(tmp, excludeState_Left));
                
                excludeState_Right = getExcludePeriods(postime,(linstate{day}{epoch}.state == 3 ));
                tmp_r = tmp(~isExcluded(tmp, excludeState_Right));
 
                posind_l = lookup(tmp_l, postime);
                posind_r = lookup(tmp_r, postime);
                
                goodpos_l = pos{day}{epoch}.data(posind_l,[2,3]);
                goodpos_r = pos{day}{epoch}.data(posind_r,[2,3]);
                %figure, plot(goodpos_l(:,1), goodpos_l(:,2),'k.');
                %figure, plot(goodpos_r(:,1), goodpos_r(:,2),'k.');
                
%                 state = find(linstate{day}{epoch}.state > 0);
%                 posinclude = intersect(ls, state);
%                 goodpos = pos{day}{epoch}.data(posinclude,:);  
                
                newpos_l{epoch} = goodpos_l;
                newpos_r{epoch} = goodpos_r;
                
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
                            tmp = tmp(~isExcluded(tmp, excludeSpeed));
                            tmp = tmp(~isExcluded(tmp, excludePos));
                            tmp = tmp(~isExcluded(tmp, excludeReward));
                            tmp = tmp(~isExcluded(tmp, excludeRipples));
                            
                            tmp_l = tmp(~isExcluded(tmp, excludeState_Left));
                            tmp_r = tmp(~isExcluded(tmp, excludeState_Right));
                            
                            
                            spikeind_l = lookup(tmp_l, spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1));
                            %numspikes_l = length(spikeind_l);
                            spikepos_l = spikes{day}{epoch}{cell(1)}{cell(2)}.data(spikeind_l,[2,3]);
                            
                            spikeind_r = lookup(tmp_r, spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1));
                            %numspikes_r = length(spikeind_r);
                            spikepos_r = spikes{day}{epoch}{cell(1)}{cell(2)}.data(spikeind_r,[2,3]);
                            
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
                minx = floor(min(allpos_l(:,1))) - 1; 
                maxx = ceil(max(allpos_l(:,1))) + 1;
                binx = (minx:binsize:maxx);
                miny = floor(min(allpos_l(:,2))) - 1;
                maxy = ceil(max(allpos_l(:,2))) + 1;
                biny = (miny:binsize:maxy);
                
                [occupancy, xticks, yticks] = hist2(allpos_l(:,1), allpos_l(:,2), binx, biny);
                
                 if ~isempty(allspikes_l)
                    %1) Calculate occupancy-normalized spikerate
                    [s, BX, BY] = hist2(allspikes_l(:,1), allspikes_l(:,2), binx, biny);
                    
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
                 
                     bad = (spikerate > 100);
                    spikerate(bad) = 0;
                    
                    occupancy = timestep*occupancy;
                    z = find(occupancy <= threshocc);
                    spikerate = spikerate(~spikerate(z));
                    spikerate_l = mean(spikerate);
                    
                 else
                     spikerate_l = 0;
                 end
                 
                 
                 %RIGHT
                minx = floor(min(allpos_r(:,1))) - 1; 
                maxx = ceil(max(allpos_r(:,1))) + 1;
                binx = (minx:binsize:maxx);
                miny = floor(min(allpos_r(:,2))) - 1;
                maxy = ceil(max(allpos_r(:,2))) + 1;
                biny = (miny:binsize:maxy);
                
                [occupancy, xticks, yticks] = hist2(allpos_r(:,1), allpos_r(:,2), binx, biny);
                
                
                if ~isempty(allspikes_r)
                    %1) Calculate occupancy-normalized spikerate
                    [s, BX, BY] = hist2(allspikes_r(:,1), allspikes_r(:,2), binx, biny);
                    
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
                    
                     bad = (spikerate > 100);
                    spikerate(bad) = 0;
                    
                    occupancy = timestep*occupancy;
                    z = find(occupancy <= threshocc);
                    spikerate = spikerate(~spikerate(z));
                    spikerate_r = mean(spikerate);
                else 
                    spikerate_r = 0;
                end
                
                    trajcomp = (spikerate_l - spikerate_r)/(spikerate_l + spikerate_r);
                    
                    SIs = [SIs; si];
                    trajcomps = [trajcomps; trajcomp];

            end
            
        end
        
    end
    
    bad = find(isnan(trajcomps));
    if ~isempty(bad)
        trajcomps(bad) = [];
        SIs(bad) = [];
    end
    
%     SIs = abs(SIs);
%     trajcomps = abs(trajcomps);
    
    figure,
    plot(SIs, trajcomps, 'k.');
    hold on
    fit = polyfit(SIs,trajcomps,1);
    plot([min(SIs), max(SIs)], polyval(fit,[min(SIs), max(SIs)] ))

    
    [CC,p] = corrcoef(SIs,trajcomps);
    R = CC(1,2)
    p = p(1,2)
end