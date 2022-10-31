       
        
        %plots the phase of all spikes from each region. can specify only phase
%locked cells, or all cells
%pulls from betaphaselock files, created with cs_phaseLocking_v3.m

%% Params

[topDir, figDir]= cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};
cellregions = {'CA1','PFC'};
betaregions = {'CA1','PFC','OB'};
plOnly = 1;
numbins = 12;

%%
for br = 1:length(betaregions)
    betaregion = betaregions{br};
    
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
            
            files = dir([animDir,'PhaseLocking\',animal,'betaphaselock_',region,'-',betaregion,'*']);
            
            for d = 1:length(files)
                load(files(d).name)
                eval(['beta_phaselock = beta_phaselock',region,';']);
                
                if plOnly == 1
                    cellfilter = '(~isempty($sph) & ($prayl < 0.05))'; %only PL cells
                else
                    cellfilter = '(~isempty($sph))'; %all cells
                end
                
                cells = evaluatefilter(beta_phaselock,cellfilter);
                
                for c = 1:size(cells,1)
                    ind = cells(c,:);
                    
                    if length(beta_phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph) > 1
                        prefphase = beta_phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.betahat;
                        allprefphase = [allprefphase; prefphase];
                    end
                end
                
            end
            
        end
        
        
        
        
        %% Get bins
        
        %allprefphase = deg2rad(allprefphase);
%         bins = -pi:(2*pi/numbins):pi;
%         [count, edges] = histcounts(allprefphase, bins);
%         
%         %pct = (count./sum(count))*100;
%         newcount = [count,count];
%         
%         bins2 = [bins - pi, bins(2:end) + pi];
%         
%         bincenters = bins2(2:end) - (bins2(2)-bins2(1))/2;
%         newbins = bincenters;
%         
        
        
        
        %% Plot

        polaraxes, hold on
        polarhistogram(allprefphase,20)
        prayl = circ_rtest(allprefphase);
        figtitle = ['Population ',region,'-',betaregion];
        title(figtitle)
        if prayl < 0.05
            pause
        end
        
%         figure, hold on
%         out = bar(newbins, newcount,1);
%         set(out, 'EdgeColor',color);
%         set(out, 'FaceColor',color);
%         plot([0 0], [0 max(newcount+2)], 'k--');
%         axis([-2*pi, 2*pi, 0 max(newcount+2)])
%         xticks([-2*pi -pi 0 pi 2*pi])
%         xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
%         xlabel('Beta Phase')
%         ylabel('Number of cells')
        
        
        
        %% Stats
%         prayl = circ_rtest(allprefphase);
%         if prayl < 0.05
%             [betahat, kappa] = circ_vmpar(allprefphase);
%             alpha = linspace(-pi, pi, numbins)';
%             [pdf] = circ_vmpdf(alpha,betahat,kappa);
%             
%             binnum = lookup(betahat,alpha);
%             pdf = pdf.*(count(binnum)/pdf(binnum));
%             newpdf = [pdf; pdf];
%             newalpha = [alpha - pi; alpha + pi];
%             %newalpha = [alpha;alpha];
%             plot(newalpha,newpdf+1,'k','LineWidth',3,'Color','k');
%             % text(1, max(newcount+1), 'p < 0.05');
%         end
        
        %continue
        %% Save Figure
        figtitle = ['Population ',region,'-',betaregion];
        title(figtitle)
        newfigfile = [figDir,'PhaseLocking\',figtitle];
        %saveas(gcf,newfigfile,'fig');
        print('-dpdf', newfigfile);
        print('-djpeg', newfigfile);
        close
    end
end