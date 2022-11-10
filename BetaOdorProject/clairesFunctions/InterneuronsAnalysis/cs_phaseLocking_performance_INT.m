%cs_phaseLocking_performance_v2
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

%Uses rank-sum to compare rather than bootstrap

clear
%close all
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};


cellregions = {'CA1','PFC'};
cellcolors=[rgbcolormap('LightSalmon'); rgbcolormap('DarkTurquoise')];

eegregions = {'CA1','PFC','OB'};

band='beta';
bandcolor=rgbcolormap('MidNightBlue');
%bandcolor=rgbcolormap('Magenta');


for cr = 1:length(cellregions)
    cellregion = cellregions{cr};
    for er=1:length(eegregions)
        
        eegregion = eegregions{er};
        zrayl = [];
        kappa = [];
        meanVecL=[];
        %allcells = [];
        for a = 1:length(animals)
            animal = animals{a};
            animDir =([topDir, animal,'Expt\',animal, '_direct\']);
            try
                load([animDir,'PhaseLocking\Interneurons\',animal,'phaselock_',band,'_',cellregion,'-',eegregion]);
            catch
                continue
            end
            phaselock_correct = phaselock;
            try
                load([animDir,'PhaseLocking\Interneurons\Incorrect\',animal,'phaselock_',band,'_incorrect_',cellregion,'-',eegregion]);
                phaselock_incorrect = phaselock;
            catch
                continue
            end
            clear phaselock;
            
            
            %get phase locked cells. use cells that are phase locked on
            %correct trials.
            cellfilt = '~isempty($sph) && $prayl < 0.05';
            
            cells = evaluatefilter(phaselock_correct,cellfilt);
            
            if isempty(cells)
                continue
            end
            
            cells = unique(cells(:,[1 3 4]),'rows');
            
            
            
            for c = 1:size(cells,1)
                ind = cells(c,:);
                %allcells = [allcells;a,ind];
                %
                k_i = phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.kappa;
                k_c = mean(phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.kappa_dist);
                
                z_i = phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.zrayl;
                z_c = mean(phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.zrayl_dist);
                
                mvl_i=phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.mvl;
                mvl_c=mean(phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.mvl_dist);
                
                
                kappa = [kappa;k_c,k_i];
                %             kappa_dist = [kappa_dist, kdist];
                
                zrayl = [zrayl;z_c,z_i];
                %zrayl_dist = [zrayl_dist, zdist];
                
                meanVecL = [meanVecL; mvl_c, mvl_i];
                
                
                %             if isempty(kappa) || isempty(zrayl)
                %                 keyboard
                %             end
            end
            
        end
        
        %
        
        %{
    % Zrayl figure
    figure
    hold on
    
    p = signrank(zrayl(:,1),zrayl(:,2));

    mn_c = mean(zrayl(:,1));
        mn_i = mean(zrayl(:,2));
        err_c = stderr(zrayl(:,1));
        err_i = stderr(zrayl(:,2));
       
        plot([r,r],[mn_c+err_c,mn_c-err_c],'k-');
        plot(r, mn_c,'ko');
        
        plot([r+0.25,r+0.25],[mn_i+err_i,mn_i-err_i],'r-');
        plot(r+0.25, mn_i,'ro');
        
        
        text(r,mn_i,['p = ',num2str(p)]);
        
        
         % Kappa figure
    figure
    hold on
    
    p = signrank(kappa(:,1),kappa(:,2));

    mn_c = mean(kappa(:,1));
        mn_i = mean(kappa(:,2));
        err_c = stderr(kappa(:,1));
        err_i = stderr(kappa(:,2));
       
        plot([r,r],[mn_c+err_c,mn_c-err_c],'k-');
        plot(r, mn_c,'ko');
        
        plot([r+0.25,r+0.25],[mn_i+err_i,mn_i-err_i],'r-');
        plot(r+0.25, mn_i,'ro');
        
        
        text(r,mn_i,['p = ',num2str(p)]);
        %}
        figure('Position',[387+360*(cr-1) 880-460*(er-1) 360 460]);
        subplot(2,1,1);
        figtitle = sprintf('Int %s, LFP %s band in %s MVL',cellregion,band,eegregion);
        z_std = std(meanVecL,1);
        z_mean = mean(meanVecL,1);
        plot(repmat([1 2],length(meanVecL),1)',meanVecL','-','color',...
            bandcolor,'LineWidth',2);
        hold on;
        if size(meanVecL,1)>=1
            
            %errorbar(z_mean,z_std,'k.')
            %axis([0 3 0 max(meanVecL(:))*1.1])
            xlim([0 3]); ylim([0 max(meanVecL(:))*1.1]);
            xticks([1 2]);
            switch band
                case 'beta'
                    ylim([0 .4]);
                case 'resp'
                     ylim([0 .75]);
            end
            xticklabels({'Correct','Incorrect'})
            [~,p] = ttest(meanVecL(:,1),meanVecL(:,2));
            if p<.05
                text(2.1, 3*(max(z_std(:))),{'ttest', ['p = ',num2str(p)]},'Color','red');
            else
                text(2.1, 3*(max(z_std(:))),{'ttest', ['p = ',num2str(p)]});
            end
            title(figtitle); box off;
            subplot(2,1,2);
            %hg=histogram(diff(meanVecL,1,2),-.15:.01:.15);
            [hc,edges]=histcounts(-diff(meanVecL,1,2),-.15:.03:.15);
            centers=mean([edges(2:end); edges(1:end-1)]);
            barh(centers,hc,1,'LineStyle','none','FaceColor',cellcolors(cr,:)); 
            hold on; plot([0 max(hc)],[0 0],'k-'); xlim([0 max(hc)]);
            %title(sprintf('ranksum p= %.2e',ranksum(meanVecL(:,1),meanVecL(:,2))));
            title(sprintf('signrank p= %.2e',signrank(meanVecL(:,1),meanVecL(:,2))));
            xlabel('Cmvl-Imvl'); ylabel('rel Proportion');
            box off; 
        end

    end
end


%{
figure(1)
title('zrayl')
xlim([0 5])
xticks([0:size(regions,1)]);
xticklabels(xlabelstr)
%axis([0 r+1 -1 9]);
% if strcmp(freq,'beta')
%     ylim([0 3])
% end
% if strcmp(freq,'resp')
%     ylim([0 30])
% end
figfile = [figDir,'Interneurons\',band,'_INT_plPerf_zrayl'];
print('-dpdf', figfile);
print('-djpeg', figfile);


figure(2)
title('kappa')
xlim([0 5])
xticks([0:size(regions,1)]);
xticklabels(xlabelstr)
%axis([0 r+1 -1 9]);
% if strcmp(freq,'beta')
%     ylim([0 3])
% end
% if strcmp(freq,'resp')
%     ylim([0 30])
% end
figfile = [figDir,'Interneurons\',band,'_INT_plPerf_kappa'];
print('-dpdf', figfile);
print('-djpeg', figfile);
%}
