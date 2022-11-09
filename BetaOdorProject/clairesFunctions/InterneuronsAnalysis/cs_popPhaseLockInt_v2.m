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
cellcolors=[rgbcolormap('LightSalmon'); rgbcolormap('DarkTurquoise')];



plOnly = 0;
numbins = 20;
freq = 'resp';
saveout=0; histbar=0;
figure;
%%
for er = 1:length(eegregions)
    eegregion = eegregions{er};
    
    for cr = 1:length(cellregions)
        region = cellregions{cr};
        
        % pull npints (unnecessary because we use the npints for the
        % calculation
        %load([topDir,'AnalysesAcrossAnimals\npInt_',region])
        
        
        
        allsph = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir = [topDir,animal,'Expt\',animal,'_direct\'];
            
            try
                files = dir([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',region,'-',eegregion,'2*']);
                % this loads the data with the CURRENT TASK RESPONSIVE SET
                load([animDir,'PhaseLocking\',files(1).name])
            catch
                fprintf('Couldnt find file for %s \n', animal)
                continue
            end
            
            if plOnly == 1
                cellfilter = '(~isempty($sph) & ($prayl < 0.05))'; %only PL cells
            else
                cellfilter = '(~isempty($sph))'; %all cells
            end
            
            
            cells = evaluatefilter(phaselock,cellfilter);
            %cells = cells(:,[1 3 4]);
            
            for c = 1:size(cells,1)
                ind = cells(c,:);
                
                if length(phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.sph) > 1
                    prefphase = phaselock{ind(1)}{ind(2)}{ind(3)}{ind(4)}.prefdir;
                    allsph = [allsph; prefphase];
                end
            end
            
        end
        
        
        %Prefphase.([region,'_',eegregion]) = allprefphase;
        
        
        % Get bins
        if histbar==1
            %allprefphase = deg2rad(allprefphase);
            bins = -pi:(2*pi/numbins):pi;
            [count, edges] = histcounts(allsph, bins);
            
            %pct = (count./sum(count))*100;
            newcount = [count,count];
            bins2 = [bins - pi, bins(2:end) + pi];
            bincenters = bins2(2:end) - (bins2(2)-bins2(1))/2;
            newbins = bincenters;
            % Plot
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
        else
            subplot(2,3,er+(cr-1)*3);
            ph=polarhistogram(allsph,numbins,'Normalization','pdf');
            set(ph,'LineStyle','none','FaceColor',cellcolors(cr,:));
            hold on;
            [rho]=circ_mean(allsph);
            [mvl]=circ_r(allsph);
            [pval, z]=circ_rtest(allsph);
            PolarArrow(rho,mvl,[],cellcolors(cr,:).*.6);
            title(sprintf('Ints %s, %s %s \n n=%d mvl=%.2f p=%.2e',cellregions{cr},...
                freq, eegregions{er},length(allsph),z, pval));
            set(gca,'RTickLabel',{},'ThetaTick',[0 90 180 270]); %,'ThetaTickLabel',{'0','\pi/2','\pi','3\pi/2'});
        end
        
        if histbar==1
            %% Stats
            [prayl,zrayl] = circ_rtest(allsph);
            nl = newline;
            text(1, max(newcount+1), ['p = ',num2str(prayl), nl, 'z = ',num2str(zrayl)]);
            if prayl < 0.05
                [hat, kappa] = circ_vmpar(allsph);
                alpha = linspace(-pi, pi, numbins)';
                [pdf] = circ_vmpdf(alpha,hat,kappa);
                
                binnum = lookup(hat,alpha);
                pdf = pdf.*(count(binnum)/pdf(binnum));
                newpdf = [pdf; pdf];
                newalpha = [alpha - pi; alpha + pi];
                %newalpha = [alpha;alpha];
                plot(newalpha,newpdf+1,'k','LineWidth',3,'Color','k');
                
            end
        end
        if saveout==1
            
            %% Save Figure
            figtitle = ['Int_',freq,' ',region,'-',eegregion];
            title(figtitle)
            newfigfile = [figDir,'PhaseLocking\PopulationCells\',figtitle];
            %saveas(gcf,newfigfile,'fig');
            print('-dpdf', newfigfile);
            print('-djpeg', newfigfile);
            %close
        end
        
    end
    
end
sgtitle(sprintf('PLonly = %d', plOnly));
%save([topDir,'AnalysesAcrossAnimals\prefPhaseInt_',freq],'Prefphase');


