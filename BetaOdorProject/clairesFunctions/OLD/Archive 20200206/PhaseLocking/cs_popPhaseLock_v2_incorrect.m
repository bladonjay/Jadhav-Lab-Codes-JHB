%plots the phase of all spikes from each region. can specify only phase
%locked cells, or all cells
%pulls from betaphaselock files, created with cs_phaseLocking_v3.m

%% Params
clear
[topDir, figDir]=cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
eegregions = {'CA1','PFC','OB'};
plOnly = 1;
numbins = 20;
freq = 'beta';

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
        allsph_test = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir,animal,'Expt\',animal,'_direct\'];
            
            files = dir([animDir,'PhaseLocking\Incorrect\',animal,freq,'phaselock_',region,'-',eegregion,'*']);
            
            for d = 1:length(files)
                load([animDir,'PhaseLocking\Incorrect\',files(d).name])
                eval(['phaselock = ',freq,'_phaselock',region,';']);
                
                %get cells that are phase locked on correct trials
                newfname = erase(files(d).name,'incorrect_');
                try
                load([animDir,'PhaseLocking\',newfname])
                catch
                    continue
                end
                eval(['phaselock_c = ',freq,'_phaselock',region,';']);
                
                %if plOnly == 1
                    cellfilter = '(~isempty($sph) & ($prayl < 0.05))'; %only PL cells
                %else
                  %  cellfilter = '(~isempty($sph))'; %all cells
                %end
                
                cells = evaluatefilter(phaselock_c,cellfilter);
                
                cellfilter = '(~isempty($sph))';
                cells_test = evaluatefilter(phaselock,cellfilter);
                
                for c = 1:size(cells,1)
                    ind = cells(c,:);
                    
                    sph = double(phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph);
                    allsph = [allsph; sph];
                end
                
%                 for c = 1:size(cells,1)
%                     ind = cells(c,:);
%                     
%                     sph_test = double(phaselock_c{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph);
%                     allsph_test = [allsph; sph_test];
%                 end
                
                
            end
            
        end
        
        if isempty(allsph)
           continue
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
        figtitle = ['PopulationSpikes ',freq,' ',region,'-',eegregion,'_incorrect'];
        title(figtitle);
        newfigfile = [figDir,'PhaseLocking\PopulationSpikes\',figtitle];
        %saveas(gcf,newfigfile,'fig');
        print('-dpdf', newfigfile);
        print('-djpeg', newfigfile);
        %close
    end
end