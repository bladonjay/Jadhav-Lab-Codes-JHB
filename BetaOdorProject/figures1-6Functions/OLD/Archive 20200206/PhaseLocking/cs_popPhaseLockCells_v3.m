%Combine all phase locked cells from CA1 and PFC

%% Params
clear
close all
[topDir, figDir]= cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
%betaregions = {'CA1','PFC','OB'};
plOnly = 1;
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
        days = cs_getRunEpochs(animDir, animal, 'odorplace');
        days = unique(days(:,1));
        
        
        files = dir([animDir,'PhaseLocking\',animal,freq,'phaselock_',region,'-','*']);
        names = {files.name};
        
        %find all cells that are phase locked
        for d = days'
            daystr = getTwoDigitNumber(d);
            dayfiles = {files(find(contains(names,daystr))).name};
            
            allcells = [];
            for j = 1:length(dayfiles)
                dayfile = dayfiles{j};
%                 if contains(dayfile, 'OB')
%                     continue
%                 end
                load([animDir,'PhaseLocking\',dayfile]);
                eval(['phaselock = ',freq,'_phaselock',region,';']);
                
                if plOnly == 1
                    cellfilter = '(~isempty($sph) & ($prayl < 0.05))'; %only PL cells
                else
                    cellfilter = '(~isempty($sph) & ~isnan($prayl))'; %all cells
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
                for f = 1:length(dayfiles)
                    dayfile = dayfiles{f};
                    load([animDir,'PhaseLocking\',dayfile]);
                    eval(['phaselock = ',freq,'_phaselock',region,';']);
                    try
                        if length(phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph) > 10 && ...
                                (phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.prayl <0.05)
                            eval(['prefphase = phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.',freq,'hat;'])
                            cellphases = [cellphases; prefphase];
                        end
                    catch
                        continue
                    end
                    
                    
                end
                
                cellprefphase = mean(cellphases);
                allprefphase = [allprefphase;cellprefphase];
                
            end
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
    prayl = circ_rtest(allprefphase);
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
    figtitle = ['Cells_',freq,'_',region,'-ALL'];
    title(figtitle)
    newfigfile = [figDir,'PhaseLocking\PopulationCells\',figtitle];
    %saveas(gcf,newfigfile,'fig');
    print('-dpdf', newfigfile);
    print('-djpeg', newfigfile);
    close
    Prefphase.(region) = allprefphase;
end

pval = circ_wwtest([Prefphase.CA1;Prefphase.PFC], [repmat(1,length(Prefphase.CA1),1);repmat(2,length(Prefphase.PFC),1)]);
%Watson-Williams, basically circular anova

mnph = circ_median(Prefphase.CA1+pi);
[~,kappa] = circ_vmpar(Prefphase.CA1);


polarplot([mnph mnph],[0 kappa])
hold on

mnph = circ_median(Prefphase.PFC+pi);
[~,kappa] = circ_vmpar(Prefphase.PFC);
polarplot([mnph mnph],[0 kappa])

annotation('textbox', [0, 0.5, 0, 0], 'string', ['p = ', num2str(pval)])

figtitle = ['Cells_',freq,'_PhasePref_stats'];
    title(figtitle)
    newfigfile = [figDir,'PhaseLocking\PopulationCells\',figtitle];
    %saveas(gcf,newfigfile,'fig');
    print('-dpdf', newfigfile);
    print('-djpeg', newfigfile);
    close
    

