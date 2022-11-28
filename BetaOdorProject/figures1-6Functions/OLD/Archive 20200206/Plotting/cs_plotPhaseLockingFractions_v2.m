function cs_plotPhaseLockingFractions_v2(dataDir,figDir,regions)

cellregions = {'CA1','PFC'};

for r = 1:length(regions)
    region = regions{r};
    load([dataDir,'betaModv2_',region,'.mat'])
    
    CA1p = [betaMod_CA1cells.prayl];
    CA1sigfract = (sum(CA1p < 0.05))/length(betaMod_CA1cells);
    
    PFCp = [betaMod_PFCcells.prayl];
    PFCsigfract = (sum(PFCp < 0.05))/length(betaMod_PFCcells);
    
    phaseLockingFractions(:,r) = [CA1sigfract; PFCsigfract];
end

save([dataDir,'phaseLockingFractions.mat'],'phaseLockingFractions');

figure,
colors = [rgb('MediumPurple'); rgb('MediumAquamarine')];
b = bar(phaseLockingFractions');
b(1).FaceColor = colors(1,:);
b(2).FaceColor = colors(2,:);

regionstr = {'CA1 beta','PFC beta','OB beta'};
set(gca,'xticklabel',regionstr)
ylabel('Fraction of Cells');
box off

set(gcf, 'Position', [50 50 1200 900]);
set(gca,'fontsize',30);




legend(cellregions)
legend('boxoff')

figtitle = 'phaseLockedFractions';
figfile = [figDir,'NicePPTFigures\',figtitle];
%saveas(gcf,figfile,'fig');
print('-bestfit','-dpdf', figfile);
print('-djpeg', figfile);

        
    
