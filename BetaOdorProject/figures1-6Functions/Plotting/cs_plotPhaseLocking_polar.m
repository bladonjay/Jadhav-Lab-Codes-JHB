close all
freq = 'beta';
[topDir,figDir] = cs_setPaths();
load([topDir,'AnalysesAcrossAnimals\prefPhaseCells_',freq]);

%get bins
bins = -pi:(2*pi/16):pi;
bins2 = bins - pi;
bincenters = bins2(2:end) - (bins2(2)-bins2(1))/2;
figure

%CA1
[count] = histcounts(Prefphase.CA1_CA1, bins,'Normalization','pdf');
smoothed = smoothdata(count,'gaussian',5);
bincenters = [bincenters,bincenters(1)];
smoothed = [smoothed,smoothed(1)];
polarplot(bincenters,smoothed,'r')

hold on
mnph = circ_median(Prefphase.CA1_CA1);
mnph = mnph-pi;
polarplot([mnph mnph],[0 0.3],'r')

%PFC
[count, edges] = histcounts(Prefphase.PFC_PFC, bins,'Normalization','pdf');
smoothed = smoothdata(count,'gaussian',5);
smoothed = [smoothed,smoothed(1)];
polarplot(bincenters,smoothed,'b')
mnph = circ_median(Prefphase.PFC_PFC);
mnph = mnph-pi;
polarplot([mnph mnph],[0 0.3],'b')

%Stats
pval = circ_wwtest([Prefphase.CA1_CA1;Prefphase.PFC_PFC], [repmat(1,length(Prefphase.CA1_CA1),1);repmat(2,length(Prefphase.PFC_PFC),1)]);

figtitle = ['Cells_',freq,'_PhasePref_polar'];
figfile = [figDir,'PhaseLocking\PopulationCells\',figtitle];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
