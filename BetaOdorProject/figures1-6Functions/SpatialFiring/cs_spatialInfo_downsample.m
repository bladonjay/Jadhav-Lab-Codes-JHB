%calculate spatial information. Compare average between odor selective
%cells, and non-selective cells.

%Use linpos/linfields- do I add the spatial information values for each
%trajectory? Or average? How does this work?

%Load selective cells and np cells, determine which cells are
%non-selective. Loop through cells, load linfields. should have occupancy
%and spikerate. Sum over epochs?

%Downsample - reduce size of nonselective cell list to equal selective
%cells (look across all animals for this. create new lists of selective and
%nonselective np cells). Get mean from new population, repeat 1000 times.
%Get dist, compare to selective (which should be the same). 

clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animalstouse = [1 2 3 4 8];

regions = {'CA1','PFC'};

[topDir, figDir] = cs_setPaths();
speedthresh = 10;
binsize = 1; % cm square
std = 3;
threshocc = 0.02; % Threshold occupancy in seconds
npExcludeDist = 30; %cm square around NP

for r = 1:length(regions)
    region = regions{r};
    
    load([topDir,'AnalysesAcrossAnimals\selectiveCells_',region]);
    load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
    
    selectivecells = selectivecells( ismember( selectivecells(:,1), animalstouse ) ,:);
    npCells = npCells( ismember( npCells(:,1), animalstouse ) ,:);
    
    nonsel = setdiff(npCells,selectivecells,'rows');
    
    %Select a random sample of nonselective cells
    nonselmeans = [];
    for I = 1:500
        samp = randsample(size(nonsel,1),size(selectivecells,1));
        nonselectivecells = nonsel(samp,:);
        disp(['Doing iteration number ', num2str(I)])
    
    spinfo_s = [];
    spinfo_n = [];
    for a = animalstouse
        animal = animals{a};
        animDir = [topDir, animal,'Expt\',animal,'_direct\'];
        
        cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
        cellfilt = ['strcmp($area,''',region,''') && strcmp($type, ''pyr'')'];
        %animCells = evaluatefilter(cellinfo,cellfilt);
        %animCells = unique(animCells(:,[1,3,4]),'rows');
        selcells = selectivecells(selectivecells(:,1)==a,2:4);
        nonselcells = nonselectivecells(nonselectivecells(:,1)==a,2:4);
        
        animCells = [selcells;nonselcells];
        animCells = sortrows(animCells,1);
        
        days = unique([selcells(:,1); nonselcells(:,1)]);
        
        for d = days'
            dayCells = animCells(animCells(:,1) == d,2:3);
            
            epochs = cs_getRunEpochs(animDir, animal, 'odorplace',d);
            
            daystr = getTwoDigitNumber(d);
            load([animDir,animal,'linfields',daystr,'.mat']);
            
            for s = 1:size(dayCells,1)
                cell = dayCells(s,:);
                spinfo = [];
                
                for ep = epochs(:,2)'
                    
                    % Linfields columns:
                    % 1) linear bin location,
                    % 2) bin occupancy (seconds),
                    % 3) bin spike count,
                    % 4) occ normailized firing per bin, and
                    % 5) smoothed occ normalized firing.
                    % 6) smoothed occupancy
                    % 7) smoothed spike count
                    
                    %Exclude cells that did not spike on track at all 
                    if length(linfields{d}{ep}) < cell(1) || length(linfields{d}{ep}{cell(1)}) < cell(2) ... 
                           || isempty(linfields{d}{ep}{cell(1)}{cell(2)})
                        continue
                    else
                        trajdata = linfields{d}{ep}{cell(1)}{cell(2)};
                        occ = cellfun(@(x) x(x(:,1)>npExcludeDist,2), trajdata([1,3]),'UniformOutput',false);
                        occ = cat(1,occ{:});
                        probocc = occ./sum(occ);
                        
                        rate = cellfun(@(x) x(x(:,1)>npExcludeDist,4), trajdata([1,3]),'UniformOutput',false);
                        rate = cat(1,rate{:});
                        normrate = rate./nanmean(rate);  %unsmoothed occ normd firing rate / mean firing rate
                        
                        if any(normrate)
                            %get rid of bins with no firing
                            incl = find(normrate);
                            normrate = normrate(incl);
                            probocc = probocc(incl);
                            
                            %get spatial info
                            spinfoperbin = probocc .* normrate .* log2(normrate); %not summed over all bins
                            
                            %sum over all bins over all traj
                        else
                            
                            spinfoperbin = NaN;
                            
                        end
                        spinfo = [spinfo, nansum(spinfoperbin)];
                        %end
                    end
                end
                mnspinfo = nanmean(spinfo);
                
                
                if ismember([d, cell], selcells, 'rows')
                    spinfo_s = [spinfo_s, mnspinfo];
                elseif ismember([d, cell], nonselcells, 'rows')
                    spinfo_n = [spinfo_n, mnspinfo];
                end
            end
        end
        
    end
    
    spinfo_s = spinfo_s(~isnan(spinfo_s));
    spinfo_n = spinfo_n(~isnan(spinfo_n));
    
    nonselmeans = [nonselmeans; mean(spinfo_n)];
    end
    
    [~,p] = ttest(nonselmeans,mean(spinfo_s))
    histogram(nonselmeans)
    %dnsmp = ceil(length(spinfo_n)/length(spinfo_s));
    %spinfo_n = downsample(sort(spinfo_n),dnsmp);
    
%     pval = cs_bootstrap(spinfo_s, spinfo_n, 10000)
%     
%     [~,p] = ttest2(spinfo_s, spinfo_n)
%     p = round(p,3);
%     stderr_s = stderr(spinfo_s);
%     stderr_n = stderr(spinfo_n);
    
    
    
    figure,
    bar([1 2],[nanmean(spinfo_s) nanmean(spinfo_n)]);
% %     hold on
% %     errorbar([nanmean(spinfo_s) nanmean(spinfo_n)], [stderr_s, stderr_n],'LineStyle','none');
    text(2,12,['p = ', num2str(p)])
    axis([0.25 2.75 0 14])
    ylabel('Spatial Information (bits/spike)');
    xticklabels({'Selective Cells','Non-selective Cells'});
    title(region)
    
    figfile = [figDir, 'PlaceFields\SpatialInfo_', region];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);
    
end