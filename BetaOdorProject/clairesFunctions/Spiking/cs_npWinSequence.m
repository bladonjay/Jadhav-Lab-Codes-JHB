%cs_npWinSequence

%plot firing rates of cells that spike in nosepoke window, ordered by
%timing of peak FR. Plot from -5 to +1.5, but only consider peaks within
%0.8 s

topDir = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

plotwin = [0 1];
%maxnp = 0.8;
binsize = 0.05;
bins = -plotwin(1):binsize:plotwin(2);

newbins = bins(1):binsize/5:bins(end-1);

regions = {'CA1','PFC'};



for r = 1:length(regions)
    region = regions{r}; 
    %load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
    load([topDir, 'AnalysesAcrossAnimals\rasterData_np_',region]);
    
    inds = rasterData.inds;
    data = rasterData.data;
    allcells = [];
    for c = 1:size(inds,1)
        celldata = data{c};
        numtrials = length(celldata.cleft) +length(celldata.cright);
        allspikes = cell2mat([celldata.cleft;celldata.cright]);
        
        counts = histcounts(allspikes,bins);
        fr = (counts/numtrials)/binsize;
        
        normfr = fr/max(fr);
        
        smoothedfr =  smoothdata(interp1(bins(1:end-1),normfr,newbins),'gaussian',15);
        allcells = [allcells;smoothedfr];
    end
    
    allcells = allcells(1:end-1,:);
    [~,peakloc] = max(allcells,[],2);
    [~,sortedinds] = sort(peakloc);
    allcells = allcells(sortedinds,:);
    figure
    imagesc(newbins,size(allcells,1):1,allcells)
    colorbar
    ylabel('Cell Number')
    xlabel('Time (seconds)')
    
    figfile = [figDir,'Spiking\npFR_',region];
       %saveas(gcf,figfile,'fig');
       print('-dpdf', figfile);
       print('-djpeg',figfile);
end

 