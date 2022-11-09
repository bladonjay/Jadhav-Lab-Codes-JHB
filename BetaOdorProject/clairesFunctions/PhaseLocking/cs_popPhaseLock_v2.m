%plots the phase of all spikes from each region. can specify only phase
%locked cells, or all cells
%pulls from betaphaselock files, created with cs_phaseLocking_v3.m

% JHB updates- now this generates dial plots and pulls from current
% betaphaselock files (postfix2)

%% Params
[topDir, figDir]=cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
cellregions = {'CA1','PFC'};
eegregions = {'CA1','PFC','OB'};
cellcolors=[rgbcolormap('LightSalmon'); rgbcolormap('DarkTurquoise')];

plOnly = 1;
numbins = 20;
freq = 'beta';
saveout=0; histbar=0;
figure;
%%
for er = 1:length(eegregions)
    eegregion = eegregions{er};
    
    for cr = 1:length(cellregions)
        region = cellregions{cr};
        
        
        
        allsph = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir,animal,'Expt\',animal,'_direct\'];
            
            files = dir([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',region,'-',eegregion,'2*']);
            
            for d = 1:length(files)
                load([animDir,'PhaseLocking\',files(d).name])
                
                if plOnly == 1
                    cellfilter = '(~isempty($sph) & ($prayl < 0.05))'; %only PL cells
                else
                    cellfilter = '(~isempty($sph))'; %all cells
                end
                
                cells = evaluatefilter(phaselock,cellfilter);
                
                for c = 1:size(cells,1)
                    ind = cells(c,:);
                    
                    sph = circ_mean(phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph);
                    allsph = [allsph; sph];
                end
                
            end
            
        end
        
        
        
        
        %% Get bins
        if histbar==1
            bins = -pi:(2*pi/numbins):pi;
            [count, edges] = histcounts(allsph, bins);
            
            pct = (count./sum(count))*100;
            newcount = [pct,pct];
            
            % going now from -2pi to 2pi- maybe to catch two cycles?
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
        else
            subplot(2,3,er+(cr-1)*3);
            ph=polarhistogram(allsph,numbins,'Normalization','pdf');
            set(ph,'LineStyle','none','FaceColor',cellcolors(cr,:));
            hold on;
            [rho]=circ_mean(allsph);
            [mvl]=circ_r(allsph);
            [pval, z]=circ_rtest(allsph);
            PolarArrow(rho,mvl,[],cellcolors(cr,:).*.6);
            title(sprintf('Cells %s, %s %s \n n=%d mvl=%.2f p=%.2e',cellregions{cr},...
                freq, eegregions{er},length(allsph),z, pval));
            set(gca,'RTickLabel',{},'ThetaTick',[0 90 180 270]); %,'ThetaTickLabel',{'0','\pi/2','\pi','3\pi/2'});
            
        end
        
        
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
        %text(1, max(newcount+1), num2str(prayl))
        
        %end
        
        %% Save Figure
        if saveout==1
            figtitle = ['PopulationSpikes ',freq,' ',region,'-',eegregion];
            title(figtitle);
            newfigfile = [figDir,'PhaseLocking\PopulationSpikes\',figtitle];
            %saveas(gcf,newfigfile,'fig');
            print('-dpdf', newfigfile);
            print('-djpeg', newfigfile);
            close
        end
    end
end
sgtitle(sprintf('PLonly = %d', plOnly));