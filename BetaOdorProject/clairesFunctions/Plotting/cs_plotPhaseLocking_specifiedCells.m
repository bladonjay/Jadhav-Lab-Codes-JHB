function cs_plotPhaseLocking_specifiedCells(animal, day, tet, cell, cellregion, eegregion, freq, savefig)
%cs_plotPhaseLocking_specifiedCells('PFC',20,'beta')
close all
nbins = 20;
[topDir, figDir] = cs_setPaths();
binedges = -pi:(2*pi/nbins):pi;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
anind = find(contains(animals,animal));

 switch cellregion
        case 'CA1'
            color = rgb('MediumIndigo');
        case 'PFC'            
            color = rgb('DarkAquamarine');
 end

    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    daystr = getTwoDigitNumber(day);
    load([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',cellregion,'-',eegregion])
    %eval(['phaselock = phaselock',cellregion,';']);

    sph = phaselock{day}{1}{tet}{cell}.sph;
    
    
    count = histcounts(sph, binedges);
                    pct = (count./sum(count))*100;
                    newcount = [pct,pct];
                    
                    newbinedges = [binedges - pi, binedges(2:end) + pi];
                    newbins = newbinedges(2:end)-((binedges(2)-binedges(1))/2);
                    
                    
                    %% Plot
                    figure, hold on
                    out = bar(newbins, newcount,1);
                    set(out, 'EdgeColor',color);
                    set(out, 'FaceColor',color);
                    plot([0 0], [0 max(newcount+2)], 'k--');
                    axis([-2*pi, 2*pi, 0 max(newcount+2)])
                    xticks([-2*pi -pi 0 pi 2*pi])
                    xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
                    xlabel('Phase')
                    ylabel('% of spikes')
                    
                    
                    
                    %% Stats
                    zrayl = phaselock{day}{1}{tet}{cell}.zrayl;
                    prayl = phaselock{day}{1}{tet}{cell}.prayl;
                    nspikes = phaselock{day}{1}{tet}{cell}.Nspikes;
                    
                    nl = newline;
                    text(0, 1, ['z = ',num2str(zrayl),nl,'p = ',num2str(prayl),nl,'n = ',num2str(nspikes)]);
    
        
    
    %% Save Figure
    if savefig == 1
                %figtitle = [region,'toCA1beta_population'];
                figtitle = [freq,' PL_',cellregion,'to',eegregion,'_',num2str(anind),'_',num2str(day),'_',num2str(tet),'_',num2str(cell)];
                    newfigfile = [figDir,'PhaseLocking\',figtitle];
                    %saveas(gcf,newfigfile,'fig');
                    print('-dpdf', newfigfile);
                    print('-djpeg', newfigfile);
    

    end
      
    
