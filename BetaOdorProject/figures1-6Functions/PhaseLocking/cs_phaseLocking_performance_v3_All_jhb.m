%cs_phaseLocking_performance_v3

%
% Jay edited to work with cs_phaseLocking current version
%
%
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

%Does this for all cells together, not split up into different eeg regions

clear
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
saveout=0;

%do all regions separately?
cellregions = {'CA1','PFC'};
cellcolors=[rgbcolormap('LightSalmon'); rgbcolormap('DarkTurquoise')];
eegregions = {'CA1','PFC','OB'};

band='beta';
bandcolor=rgbcolormap('MidNightBlue'); % this is for beta
%bandcolor=rgbcolormap('Magenta'); % this is for rr
for cr = 1:length(cellregions)
    
    cellregion = cellregions{cr};
    for er=1:length(eegregions)
        eegregion=eegregions{er};
        zrayl = [];
        kappa = [];
        meanVecL = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir =([topDir, animal,'Expt\',animal, '_direct\']);
            days = cs_getRunEpochs(animDir, animal, 'odorplace');
            days = unique(days(:,1));
            
            % phaselocking correct data (pull all NP cells here for all
            % cells
            plfiles_c = dir([animDir,'PhaseLocking\',animal,'phaselock_', band, '_', cellregion, '-', eegregion, '2*']);
            names_c = {plfiles_c.name};
            % phaselocking incorrect data
            plfiles_i = dir([animDir,'PhaseLocking\Incorrect\',animal,'phaselock_', band, '_',...
                'incorrect_', cellregion, '-', eegregion, '2*']);
            names_i = {plfiles_i.name};
            if length(plfiles_c)~=1 || length(plfiles_i)~=1
                continue
            end
            
            % load up the raw data
            load([animDir,'PhaseLocking\',plfiles_c.name]);
            eval('phaselock_c = phaselock;');
            clear phaselock;
            load([animDir,'PhaseLocking\Incorrect\',plfiles_i.name]);
            eval('phaselock_i = phaselock;');
            clear phaselock;
            
            %get phase locked cells. use cells that are phase locked on
            %correct trials.
            %cellfilt = '~isempty($sph) & $prayl <0.05';
            cellfilt = '~isempty($sph)';
            allcells = evaluatefilter(phaselock_c,cellfilt);
            if isempty(allcells)
                continue
            end
            cells = unique(allcells(:,[1 3 4]),'rows');
            %For each cell, get avg kappa and zrayl.
            for c = 1:size(cells,1)
                ind=cells(c,:);
                
                
                %get kappa and zrayl scores for comparison
                k_ep = [mean(phaselock_c{ind(1)}{1}{ind(2)}{ind(3)}.kappa_dist), phaselock_i{ind(1)}{1}{ind(2)}{ind(3)}.kappa];
                z_ep = [mean(phaselock_c{ind(1)}{1}{ind(2)}{ind(3)}.zrayl_dist), phaselock_i{ind(1)}{1}{ind(2)}{ind(3)}.zrayl];
                mvl_ep=[mean(phaselock_c{ind(1)}{1}{ind(2)}{ind(3)}.mvl_dist), phaselock_i{ind(1)}{1}{ind(2)}{ind(3)}.mvl];
                
                kappa = [kappa; k_ep];
                zrayl = [zrayl;z_ep];
                meanVecL = [meanVecL; mvl_ep];
                
                
            end
        end
        
        %%
        % Kappa figure
        %{
        figure;
        k_std = std(kappa,1);
        k_mean = mean(kappa,1);
        plot(repmat([1 2],length(kappa),1)',kappa','-');
        hold on
        errorbar(k_mean,k_std,'k.')
        figtitle = sprintf('Pyrs %s, LFP %s kappa',cellregion,eegregion);
        k_std = std(kappa,1);
        k_mean = mean(kappa,1);
        errorbar(k_mean,k_std,'k.')
        axis([0 3 0 (max(k_mean(:)) +2*(max(k_std(:))))])
        xticks([1 2])
        xticklabels({'Correct','Incorrect'})
        [~,p] = ttest(kappa(:,1),kappa(:,2));
        text(1.5, 3*(max(k_std(:))),['p = ',num2str(p)]);
        title(figtitle)
        if saveout==1
            figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
        end
        % Zrayl figure
        figure;
        figtitle = sprintf('Pyrs %s, LFP %s zrayl',cellregion,eegregion);
        z_std = std(zrayl,1);
        z_mean = mean(zrayl,1);
        plot(repmat([1 2],length(zrayl),1)',zrayl','-');
        hold on;
        errorbar(z_mean,z_std,'k.')
        axis([0 3 0 (max(z_mean(:)) +2*(max(z_std(:))))])
        xticks([1 2])
        xticklabels({'Correct','Incorrect'})
        [~,p] = ttest(zrayl(:,1),zrayl(:,2));
        text(1.5, 3*(max(z_std(:))),['p = ',num2str(p)]);
        title(figtitle)
        if saveout==1
            figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
        end
        %}
        % MVL figure
        figure('Position',[1474+360*(cr-1) 880-460*(er-1) 360 460]);

        subplot(2,1,1);
        figtitle = sprintf('Pyrs %s, LFP %s band in %s MVL',cellregion,band,eegregion);
        z_std = std(meanVecL,1);
        z_mean = mean(meanVecL,1);
        plot(repmat([1 2],length(meanVecL),1)',meanVecL','-','color',...
            bandcolor,'LineWidth',2);
        hold on;
        %errorbar(z_mean,z_std,'k.')
        %axis([0 3 0 (max(z_mean(:)) +2*(max(z_std(:))))])
        xlim([0 3]); ylim([0 .75]);
        xticks([1 2])
        xticklabels({'Correct','Incorrect'})
        [~,p] = ttest(meanVecL(:,1),meanVecL(:,2));
        if p<.05
           text(2.1, 3*(max(z_std(:))),{'ttest', ['p = ',num2str(p)]},'Color','red');
        else
           text(2.1, 3*(max(z_std(:))),{'ttest', ['p = ',num2str(p)]});
        end
        title(figtitle); box off
        subplot(2,1,2);
        % col1 minus col2
        [hc,edges]=histcounts(-diff(meanVecL,1,2),-.2:.03:.2);
        centers=mean([edges(2:end); edges(1:end-1)]);
        %title(sprintf('ranksum p= %.2e',ranksum(meanVecL(:,1),meanVecL(:,2))));
        barh(centers,hc,1,'LineStyle','none','FaceColor',cellcolors(cr,:)); 
        hold on; plot([0 max(hc)],[0 0],'k-'); set(gca,'XLim',[0 max(hc)]);
        title(sprintf('signrank p= %.2e',signrank(meanVecL(:,1)-meanVecL(:,2))));
        xlabel('Cmvl-Imvl'); ylabel('rel Proportion');
        box off;
        if saveout==1
            figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
        end
        
    end
end