%plots the phase of all spikes from each region. can specify only phase
%locked cells, or all cells
%most updated version as of 10/6/2020
%% Params
clear
close all
[topDir, figDir]= cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
eegregions = {'CA1','PFC','OB'};
plOnly = 0;
numbins = 12;
freq = 'resp';
%%
for br = 1:length(eegregions)
    eegregion = eegregions{br};
    
    for cr = 1:length(cellregions)
        region = cellregions{cr};
        
        load([topDir,'AnalysesAcrossAnimals\npCells_',region])
        
        switch region
            case 'CA1'
                color = rgb('MediumPurple');
            case 'PFC'
                color = rgb('DarkAquamarine');
        end
        
        allprefphase = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir,animal,'Expt\',animal,'_direct\'];
            
            load([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',region,'-',eegregion])
            
            if plOnly == 1
                cellfilter = '(~isempty($sph) & ($prayl < 0.05))'; %only PL cells
            else
                cellfilter = '(~isempty($sph))'; %all cells
            end
            
            %cells = npCells(npCells(:,1) == a,2:4);
            
            cells = evaluatefilter(phaselock,cellfilter);
            cells = cells(:,[1 3 4]);
            
            for c = 1:size(cells,1)
                ind = cells(c,:);

                if length(phaselock{ind(1)}{1}{ind(2)}{ind(3)}.sph) > 1
                    prefphase = phaselock{ind(1)}{1}{ind(2)}{ind(3)}.prefdir;
                    allprefphase = [allprefphase; prefphase];
                end
            end
            
        end
        
        
        Prefphase.([region,'_',eegregion]) = allprefphase;
        
        
        %% Get bins
        
        %allprefphase = deg2rad(allprefphase);
        bins = -pi:(2*pi/numbins):pi;
        [count, edges] = histcounts(allprefphase, bins);
        
        %pct = (count./sum(count))*100;
        newcount = [count,count];
        
        bins2 = [bins - pi, bins(2:end) + pi];
        
        bincenters = bins2(2:end) - (bins2(2)-bins2(1))/2;
        newbins = bincenters;
        
        
        
        
        %% Plot
        
        figure, hold on
        out = bar(newbins, newcount,1);
        set(out, 'EdgeColor',color);
        set(out, 'FaceColor',color);
        plot([0 0], [0 max(newcount+2)], 'k--');
        axis([-2*pi, 2*pi, 0 max(newcount+2)])
        xticks([-2*pi -pi 0 pi 2*pi])
        xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
        xlabel([freq,' phase'])
        ylabel('Number of cells')

        
        
        %% Stats
        [prayl,zrayl] = circ_rtest(allprefphase);
        nl = newline;
        text(1, max(newcount+1), ['p = ',num2str(prayl), nl, 'z = ',num2str(zrayl)]);
        if prayl < 0.05
            [hat, kappa] = circ_vmpar(allprefphase);
            alpha = linspace(-pi, pi, numbins)';
            [pdf] = circ_vmpdf(alpha,hat,kappa);
            
            binnum = lookup(hat,alpha);
            pdf = pdf.*(count(binnum)/pdf(binnum));
            newpdf = [pdf; pdf];
            newalpha = [alpha - pi; alpha + pi];
            %newalpha = [alpha;alpha];
            plot(newalpha,newpdf+1,'k','LineWidth',3,'Color','k');
            
        end
        
        
        %% Save Figure
        figtitle = ['Cells_',freq,' ',region,'-',eegregion];
        title(figtitle)
        newfigfile = [figDir,'PhaseLocking\PopulationCells\',figtitle];
        %saveas(gcf,newfigfile,'fig');
        print('-dpdf', newfigfile);
        print('-djpeg', newfigfile);
        %close
        
        
    end
    
end
 save([topDir,'AnalysesAcrossAnimals\prefPhaseCells_',freq],'Prefphase');
 
 
