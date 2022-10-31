%function cs_plotPhaseLocking_specifiedCells(region,nbins)
clear
close all
nbins = 20;
regions = {'CA1','PFC'};
eegregions = {'CA1','PFC'};
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
binedges = -pi:(2*pi/nbins):pi;
freq = 'beta';

for cr = 1:length(regions)
    region = regions{cr};
    switch region
        case 'CA1'
            color = rgb('MediumIndigo');
        case 'PFC'
            color = rgb('DarkAquamarine');
    end
    
    
    load([topDir, 'AnalysesAcrossAnimals\PhaseLocking\plCells_',freq,'_',region,'-',region,'.mat']);
    
    aninds = unique(plcells(:,1));
    
    
    for a = 1:length(aninds)
        anind = aninds(a);
        animal = animals{anind};
        animDir = [topDir, animal,'Expt\',animal,'_direct\'];
        
        ancells = plcells(plcells(:,1) == anind,2:4);
        
        
        load([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',region,'-',region,'.mat'])
        
        
        for c = 1:size(ancells,1)
            cell = ancells(c,:);
            
            %cellfilter = '(~isempty($sph) & ($prayl < 0.05) &($Nspikes > 10))'; %only PL cells
            %testcells = evaluatefilter(phaselock{day}{epoch},cellfilter);
            sph = phaselock{cell(1)}{1}{cell(2)}{cell(3)}.sph;
            
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
            
            
            zrayl = phaselock{cell(1)}{1}{cell(2)}{cell(3)}.zrayl;
            prayl = phaselock{cell(1)}{1}{cell(2)}{cell(3)}.prayl;
            nl = newline;
            text(0, 1, ['z = ',num2str(zrayl),nl,'p = ',num2str(prayl)]);
            
            
            %title([animal, ' ', num2str
            %% Save Figure
            figtitle = [freq,' PL_',region,'to',region,'_',num2str(anind),'_',num2str(cell(1)),'_',num2str(cell(2)),'_',num2str(cell(3))];
            newfigfile = [figDir,'PhaseLocking\',figtitle];
            %saveas(gcf,newfigfile,'fig');
            print('-dpdf', newfigfile);
            print('-djpeg', newfigfile);
            
            close
        end
        
    end
end

