%cs_npWinSequence

%plot firing rates of cells that spike in nosepoke window, ordered by
%timing of peak FR. Plot from -5 to +1.5, but only consider peaks within
%0.8 s

topDir = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

plotwin = [0.5 1.5];
%maxnp = 0.8;
binsize = 0.05;
plotbins = -plotwin(1):binsize:plotwin(2);
trigbin = plotwin(1)/binsize;


regions = {'CA1'};
trialtypes = {'alltrials','lefttrials','righttrials'};

for r = 1:length(regions)
    region = regions{r}; 
    load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
    animNums = unique(npCells(:,1));
    allpeaks = [];
    allplotfr = [];
    for a = animNums'
        animal = animals{a};
        animDir = [topDir, animal, 'Expt/',animal,'_direct/'];
        animCells = npCells(npCells(:,1) == a,2:4);
        days = unique(animCells(:,1));
        epochs = cs_getRunEpochs(animDir, animal, 'odorplace');
        cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
        
        for day = days'
            
            disp(['Animal ',animal,' day ',num2str(day)])
            daycells = animCells(animCells(:,1)==day, 2:3);
            
            pyrcellfilt = 'strcmp($type,''pyr'')'; %pyramidal cells only
            pyrcells = evaluatefilter(cellinfo{day},pyrcellfilt);
            pyrcells = unique(pyrcells(:,2:3),'rows');
            daycells = intersect(daycells,pyrcells,'rows');
            
            npWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
            spikes = loaddatastruct(animDir, animal, 'spikes',day);
            dayeps = epochs(epochs(:,1)==day,2);
            
            for c = 1:size(daycells,1)
                cell = daycells(c,:);

                allnpwinspikes = [];
                allplotspikes = [];
                totaltrials = 0;
                for epoch = dayeps'
  
                    wins = npWindow{day}{epoch};
                    if ~isempty(spikes{day}{epoch}{cell(1)}{cell(2)})
                    sp = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                    totaltrials = totaltrials + size(wins,1);
                    else
                        continue
                    end
                    %calculate peak np spiking time from actual NP window
                    %spikes. Then, collect spike times that fall within
                    %plot window
                    for w = 1:size(wins,1)
                        window = wins(w,:);
                        winsp = sp(isExcluded(sp, window)); 
                        %isExlcuded: 1 signifies times that fall in between start and
                        %end times of excludePeriods
                        winsp = winsp - window(1);
                        
                        allnpwinspikes = [allnpwinspikes;winsp];
                        
                        plotwindow = [window(1)-plotwin(1), window(1)+plotwin(2)];
                        plotsp = sp(isExcluded(sp, plotwindow)); 
                        plotsp = plotsp - window(1);
                        
                        allplotspikes = [allplotspikes;plotsp];
                    end
                    
                end
                counts = histcounts(allnpwinspikes,15);
                if sum(counts)< totaltrials
                    keyboard
                end
                avgcounts = counts./totaltrials;
                
                smoothedcounts = smoothdata(avgcounts,'gaussian',3);
                [~,peakind] = max(smoothedcounts);
                
                
                allpeaks = [allpeaks,peakind];
                %get average fr
                plotcounts = histcounts(allplotspikes, plotbins);
                avgplotcounts = plotcounts./totaltrials;
                plotfr = avgplotcounts./binsize;
                plotfr = smoothdata(plotfr,'gaussian',5);
                
                allplotfr = [allplotfr;plotfr];
            end
        end
    end
    
    [~,inds] = sort(allpeaks);
    allplotfr = allplotfr(inds, :);
    
    figure, 
    imagesc(allplotfr)
    colormap jet
    hold on
    plot([trigbin, trigbin], [1 size(allplotfr,1)], 'k--')
    colorbar
    
end