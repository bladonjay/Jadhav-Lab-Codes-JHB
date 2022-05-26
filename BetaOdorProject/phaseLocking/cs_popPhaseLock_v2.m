%plots the phase of all spikes from each region. can specify only phase
%locked cells, or all cells
%pulls from betaphaselock files, created with cs_phaseLocking_v3.m

%% Params
[topDir, figDir]=cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
eegregions = {'OB','CA1','PFC'};
plOnly = 1;
numbins = 20;
freq = 'resp';

%%
for er = 1:length(eegregions)
    eegregion = eegregions{er};
    
    for cr = 1:length(cellregions)
        region = cellregions{cr};
        
        switch region
            case 'CA1'
                color = rgb('MediumIndigo');
            case 'PFC'
                color = rgb('DarkAquamarine');
        end
        
        allsph = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir,animal,'Expt\',animal,'_direct\'];
            
            files = dir([animDir,'PhaseLocking\',animal,freq,'phaselock_',region,'-',eegregion,'*']);
            
            for d = 1:length(files)
                load([animDir,'PhaseLocking\',files(d).name])
                eval(['phaselock = ',freq,'_phaselock',region,';']);
                
                if plOnly == 1
                    cellfilter = '(~isempty($sph) & ($prayl < 0.05))'; %only PL cells
                else
                    cellfilter = '(~isempty($sph))'; %all cells
                end
                
                cells = evaluatefilter(phaselock,cellfilter);
                
                for c = 1:size(cells,1)
                    ind = cells(c,:);
                    
                    sph = phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph;
                    allsph = [allsph; sph];
                end
                
            end
            
        end
        
        
        
        
        %% Get bins
        
        bins = -pi:(2*pi/numbins):pi;
        [count, edges] = histcounts(allsph, bins);
        
        pct = (count./sum(count))*100;
        newcount = [pct,pct];
        
        bins2 = [bins - pi, bins(2:end) + pi];
        
        bincenters = bins2(2:end) - (bins2(2)-bins2(1))/2;
        newbins = bincenters;
        
        
        
        
%% Plot
        figure, hold on
        out = bar(newbins, newcount,1);
        set(out, 'EdgeColor',color);
        set(out, 'FaceColor',color);
        plot([0 0], [0 max(newcount+2)], 'k--');
        axis([-2*pi, 2*pi, 1 max(newcount+2)])
        xticks([-2*pi -pi 0 pi 2*pi])
        xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
        xlabel('Phase')
        ylabel('% of spikes')
        
        
        
        % Stats
        [prayl] = circ_rtest(allsph);
        %if prayl < 0.05
%             [hat, kappa] = circ_vmpar(allsph);
%             alpha = linspace(-pi, pi, numbins)';
%             [pdf] = circ_vmpdf(alpha,hat,kappa);
%             
%             binnum = lookup(hat,alpha);
%             pdf = pdf.*(pct(binnum)/pdf(binnum));
%             newpdf = [pdf; pdf];
%             newalpha = [alpha - pi; alpha + pi];
%             %newalpha = [alpha;alpha];
%             plot(newalpha,newpdf+1,'k','LineWidth',3,'Color','k');
            text(1, max(newcount+1), num2str(prayl))
            
        %end
           
%% Save Figure
        figtitle = ['PopulationSpikes ',freq,' ',region,'-',eegregion];
        title(figtitle);
        newfigfile = [figDir,'PhaseLocking\PopulationSpikes\',figtitle];
        %saveas(gcf,newfigfile,'fig');
        print('-dpdf', newfigfile);
        print('-djpeg', newfigfile);
        close
    end
end