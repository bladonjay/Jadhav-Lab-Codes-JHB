%function cs_plotPhaseLocking_specifiedCells(region,nbins)
close all
nbins = 20;
regions = {'CA1','PFC'};
eegregions = {'CA1','PFC','OB'};
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
binedges = -pi:(2*pi/nbins):pi;
freq = 'resp';

for cr = 1:length(regions)
    region = regions{cr};
    switch region
        case 'CA1'
            color = rgb('MediumIndigo');
        case 'PFC'
            color = rgb('DarkAquamarine');
    end
    
    for er = 1:length(eegregions)
        eegregion = eegregions{er};
        
        load([topDir, 'AnalysesAcrossAnimals\PhaseLocking\plCells_',freq,'_',region,'-',eegregion,'.mat']);
        
        aninds = unique(plcells(:,1));
        
        
        for a = 1:length(aninds)
            anind = aninds(a);
            animal = animals{anind};
            animDir = [topDir, animal,'Expt\',animal,'_direct\'];
            ancells = plcells(plcells(:,1) == anind,2:4);
            
            days = unique(ancells(:,1));
            %files = dir([animDir,'PhaseLocking\',animal,'betaphaselock_',region,'-',eegregion,'_0*']);
            
            for d = 1:length(days)
                day = days(d);
                cells = ancells(ancells(:,1) == day,2:3);
                daystr = getTwoDigitNumber(day);
                load([animDir,'PhaseLocking\',animal,freq,'phaselock_',region,'-',eegregion,'_',daystr,'.mat'])
                
                eval(['phaselock = ',freq,'_phaselock',region,';']);
                
                epochs = find(~cellfun(@isempty, phaselock{day}));
                %             cellfilter = '(~isempty($sph) & ($prayl < 0.05) &($Nspikes > 10))'; %only PL cells
                %             cells = evaluatefilter(beta_phaselock,cellfilter);
                %
                for c = 1:size(cells,1)
                    cell = cells(c,:);
                    allsph = [];
                    for ep = 1:length(epochs)
                        epoch = epochs(ep);
                        
                        cellfilter = '(~isempty($sph) & ($prayl < 0.05) &($Nspikes > 10))'; %only PL cells
                        testcells = evaluatefilter(phaselock{day}{epoch},cellfilter);
                        if ismember(cell, testcells, 'rows')
                            
                            sph = phaselock{day}{epoch}{cell(1)}{cell(2)}.sph;
                            %                 alpha = beta_phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.alpha;
                            %                 betahat = beta_phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.betahat;
                            %                 kappa = beta_phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.kappa;
                            
                            
                            allsph = [allsph; sph];
                        end
                    end
                    sph = allsph;
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
                    
                    
                    [betahat, kappa] = circ_vmpar(sph);
                    alpha = linspace(-pi, pi, nbins)';
                    
                    [pdf] = circ_vmpdf(alpha,betahat,kappa);
                    binnum = lookup(betahat,alpha);
                    pdf = pdf.*(pct(binnum)/pdf(binnum));
                    newpdf = [pdf; pdf];
                    newalpha = [alpha - pi; alpha + pi];
                    %newalpha = [alpha;alpha];
                    %plot(newalpha,newpdf+1,'k','LineWidth',3,'Color','k');
                    
                    
                    %title([animal, ' ', num2str
                    %% Save Figure
                    figtitle = [freq,' PL_',region,'to',eegregion,'_',num2str(anind),'_',num2str(day),'_',num2str(cell(1)),'_',num2str(cell(2))];
                    newfigfile = [figDir,'PhaseLocking\',figtitle];
                    %saveas(gcf,newfigfile,'fig');
                    print('-dpdf', newfigfile);
                    print('-djpeg', newfigfile);
                    
                    close
                end
                
            end
        end
    end
    
end


%