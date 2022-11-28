%Combine all phase locked cells from CA1 and PFC

%% Params
clear
close all
[topDir, figDir]= cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
%betaregions = {'CA1','PFC','OB'};
plOnly = 0;
numbins = 12;

freq = 'beta';

for cr = 1:length(cellregions)
    region = cellregions{cr};
    
    switch region
        case 'CA1'
            color = rgb('MediumIndigo');
        case 'PFC'
            color = rgb('DarkAquamarine');
    end
    
    allprefphase = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        files = dir([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',region,'-*']);
        names = {files.name};
        
        %find all cells that are phase locked
        files = {files(find(~contains(names,'prelearn') & ~contains(names,'postlearn'))).name};
        
        allcells = [];
        for j = 1:length(files)
            file = files{j};
            %                 if contains(dayfile, 'OB')
            %                     continue
            %                 end
            load([animDir,'PhaseLocking\',file]);
            
            
            if plOnly == 1
                cellfilter = '$Nspikes > 10 & $prayl < 0.05'; %only PL cells
            else
                cellfilter = '$Nspikes > 10'; %all cells
            end
            
            cells = evaluatefilter(phaselock,cellfilter);
            allcells = [allcells;cells];
            
        end
        cells = unique(allcells,'rows');
        
        %loop over cells to find avg perferred phase
        for c = 1:size(cells,1)
            ind = cells(c,:);
            
            %go through each beta region
            cellphases = [];
            for f = 1:length(files)
                file = files{f};
                load([animDir,'PhaseLocking\',file]);
                try
                    %if length(phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph) > 10
                   prefphase = phaselock{ind(1)}{1}{ind(3)}{ind(4)}.prefdir;
                    if isnan(prefphase)
                        continue
                    end
                    cellphases = [cellphases; prefphase];
                    %end
                catch
                    continue
                end
                
                
            end
            
            cellprefphase = nanmean(cellphases);
            if isnan(cellprefphase)
                keyboard
            end
            allprefphase = [allprefphase;cellprefphase];
            
        end
    end
    
    
    
    %end
    
    
    
    
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
    text(1, max(newcount+1), ['p = ',num2str(prayl)]);
    if prayl < 0.05
        [betahat, kappa] = circ_vmpar(allprefphase);
        alpha = linspace(-pi, pi, numbins)';
        [pdf] = circ_vmpdf(alpha,betahat,kappa);
        
        binnum = lookup(betahat,alpha);
        pdf = pdf.*(count(binnum)/pdf(binnum));
        newpdf = [pdf; pdf];
        newalpha = [alpha - pi; alpha + pi];
        %newalpha = [alpha;alpha];
        %plot(newalpha,newpdf+1,'k','LineWidth',3,'Color','k');
        
        
    end
    
    
    %% Save Figure
    figtitle = ['Cells-',freq,'-',region,'-ALL'];
    title(figtitle)
    newfigfile = [figDir,'PhaseLocking\PopulationCells\',figtitle];
    %saveas(gcf,newfigfile,'fig');
    print('-dpdf', newfigfile);
    print('-djpeg', newfigfile);
    %close
    Prefphase.(region) = allprefphase;
end

% save([topDir,'AnalysesAcrossAnimals\prefPhaseCells_',freq],'Prefphase');
% pval = circ_wwtest([Prefphase.CA1;Prefphase.PFC], [repmat(1,length(Prefphase.CA1),1);repmat(2,length(Prefphase.PFC),1)]);
% %Watson-Williams, basically circular anova
% 
% mnph = circ_median(Prefphase.CA1);
% [~,kappa] = circ_vmpar(Prefphase.CA1);
% 
% %polarhistogram(Prefphase.CA1,16,'DisplayStyle','stairs')
% 
% polarplot([mnph mnph],[0 kappa])
% hold on
% %polarhistogram(Prefphase.PFC,12,'DisplayStyle','stairs')
% mnph = circ_median(Prefphase.PFC);
% [~,kappa] = circ_vmpar(Prefphase.PFC);
% polarplot([mnph mnph],[0 kappa])
% 
% annotation('textbox', [0, 0.5, 0, 0], 'string', ['p = ', num2str(pval)])
% 
% figtitle = ['Cells-',freq,'PhasePref stats'];
% title(figtitle)
% newfigfile = [figDir,'PhaseLocking\PopulationCells\',figtitle];
% %saveas(gcf,newfigfile,'fig');
% print('-dpdf', newfigfile);
% print('-djpeg', newfigfile);
% close
% 
