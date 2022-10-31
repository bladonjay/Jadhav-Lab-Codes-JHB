[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC'};
for r = 1:length(regions)
    region = regions{r};
cells = [];
for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    epochs = cs_getRunEpochs(animDir, animal, 'odorplace');
    
    days = unique(epochs(:,1));
    for d = 1:length(days)
        day = days(d);
        daystr = getTwoDigitNumber(day);
        
        load([animDir,animal,'cellinfo.mat']);
        cellfilter = ['($numspikes > 0) & strcmp($area,''',region,''')']; 
        daycells = evaluatefilter(cellinfo{day},cellfilter);
        
        uniquecells = unique(daycells(:,[2,3]),'rows');
        numcells = size(uniquecells,1);

        cells = [cells; numcells]; %#ok<AGROW>
       
    end
    
end

histogram(cells,10)
[n,g]= histcounts(cells,10);
xlabel('Number of Cells per Day')
yticks(min(n):max(n)+1)
axis([min(g)-2, max(g)+2, 0, max(n)+1])

figtitle = [region,' Number of Cells Distribution'];
figfile = [figDir,'Plots\',figtitle];
saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);
end