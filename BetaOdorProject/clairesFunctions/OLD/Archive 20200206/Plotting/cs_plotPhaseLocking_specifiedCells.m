function cs_plotPhaseLocking_specifiedCells(region,nbins,freq)
%cs_plotPhaseLocking_specifiedCells('PFC',20,'beta')

[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
bins = -pi:(2*pi/nbins):pi;

 switch region
        case 'CA1'
            color = rgb('MediumIndigo');
        case 'PFC'            
            color = rgb('DarkAquamarine');
 end

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    files = dir([animDir,'PhaseLocking\',animal,freq,'phaselock_',region,'-*']);
    
    for d = 1:length(files)
            load([animDir,'PhaseLocking\',files(d).name])
            
            eval(['phaselock = ',freq,'_phaselock',region,';']);
            eegregion = extractBefore(extractAfter(files(d).name,[animal,freq,'phaselock_',region,'-',]),'_');
            day = length(phaselock);
            %get eeg region
           
            cellfilter = '(~isempty($sph) & ($prayl < 0.05) &($Nspikes > 10))'; %only PL cells
            cells = evaluatefilter(phaselock,cellfilter);
            if ~isempty(cells)
            cells = unique(cells(:,[1, 3, 4]),'rows');
            
            
            for c = 1:size(cells,1)
                ind = cells(c,:);
                epochs = cs_findGoodEpochs(phaselock{ind(1)},{'sph'},ind(2:3));
                allsph = [];
                for ep = epochs'
                    
                sph = phaselock{ind(1)}{ep}{ind(2)}{ind(3)}.sph;
                allsph = [allsph;sph];
                end
                
                count = histc(allsph, bins);
                
                
                pct = (count./sum(count))*100;
                newcount = [pct(1:end-1);pct];
                
                newbins = [bins - pi, bins(2:end) + pi];


                %% Plot
                figure, hold on
                out = bar(newbins, newcount,1);
                set(out, 'EdgeColor',color);
                set(out, 'FaceColor',color);
                plot([0 0], [0 max(newcount+2)], 'k--');
                axis([-2*pi, 2*pi, 1 max(newcount+2)])
                xticks([-2*pi -pi 0 pi 2*pi])
                xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})
                xlabel('Beta Phase')
                ylabel(' of spikes')
    
    
    
%% Stats 
    
                  
                 [betahat, kappa] = circ_vmpar(sph);
                 alpha = linspace(-pi, pi, nbins+1)';

                 [pdf] = circ_vmpdf(alpha,betahat,kappa);
                 binnum = lookup(betahat,alpha);
                 pdf = pdf.*(pct(binnum)/pdf(binnum));
                 newpdf = [pdf(1:end-1); pdf]; 
                 newalpha = [alpha(1:end-1) - pi; alpha + pi]; 
                 %newalpha = [alpha;alpha];
                 %plot(newalpha,newpdf+1,'k','LineWidth',3,'Color','k');
        
    
    %% Save Figure
                %figtitle = [region,'toCA1beta_population'];
                figtitle = [freq,' PL_',region,'to',eegregion,'_',num2str(a),'_',num2str(day),'_',num2str(ind(2)),'_',num2str(ind(3))];
                newfigfile = [figDir,'PhaseLocking\',figtitle];
                saveas(gcf,newfigfile,'fig');
                print('-dpdf', newfigfile);
                print('-djpeg', newfigfile);
    

            end
            end
           
    end
end
    
    
    
