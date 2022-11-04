%cs_phaseLocking_performance
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

%Uses rank-sum to compare rather than bootstrap

clear
close all
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

freq = 'beta';

%get all region combinations
% [A,B] = meshgrid(cellregions,eegregions);
% c=cat(2,A',B');
% regions = reshape(c,[],2);
regions = {'CA1','CA1';'PFC','PFC'};
xlabelstr = [];
for r = 1:size(regions,1)
    cellregion = regions{r,1};
    eegregion = regions{r,2};
    
    zrayl = [];
    kappa = [];
    kappa_dist = [];
    zrayl_dist = [];
    R_all = [];
    r_dist = [];
    allcells = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir =([topDir, animal,'Expt\',animal, '_direct\']);
        load([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',cellregion,'-',eegregion]);
        phaselock_correct = phaselock;
        try
            load([animDir,'PhaseLocking\Incorrect\',animal,'phaselock_',freq,'_incorrect_',cellregion,'-',eegregion]);
            phaselock_incorrect = phaselock;
        catch
            continue
        end
        
        
        
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
            allcells = [allcells;a,ind];
            %
%             k = phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.kappa;
%             kdist = phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.kappa_dist;
            
            z_i = phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.zrayl;
            z_c = mean(phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.zrayl_dist);
            
            
            
%             kappa = [kappa;k];
%             kappa_dist = [kappa_dist, kdist];
            
            zrayl = [zrayl;z_c,z_i];
            %zrayl_dist = [zrayl_dist, zdist];
            
            
            
            
%             if isempty(kappa) || isempty(zrayl)
%                 keyboard
%             end
        end
        
    end
   
%     kappa_dist = mean(kappa_dist,2);
%     zrayl_dist = mean(zrayl_dist,2);
%     % Kappa figure
%     figure(1)
%     hold on
%     p_upper  = sum(kappa_dist >= mean(kappa)) /length(kappa_dist);
%     p_lower = sum(kappa_dist <= mean(kappa))/length(kappa_dist);
%     p =   min([p_upper, p_lower]);
%     
%     CI = getCI(kappa_dist,95);
%     plot([r r], CI','k-')
%     hold on
%     plot(r,mean(kappa_dist),'ko')
%     plot(r+0.25,mean(kappa),'ro');
%     
%     %     k_std = std(kappa,1);
%     %     k_mean = mean(kappa,1);
%     %     errorbar([r, r+0.25], k_mean,k_std,'k.')
%     xlabelstr = [xlabelstr,{[cellregion,'-',eegregion]}];
%     %     figtitle = [cellregion,'-',eegregion,' kappa'];
%     %     k_std = std(kappa,1);
%     %     k_mean = mean(kappa,1);
%     %     errorbar(k_mean,k_std,'k.')
%     %     axis([0 3 0 (max(k_mean(:)) +2*(max(k_std(:))))])
%     %     xticks([1 2])
%     %xticklabels({'Correct','Incorrect'})
%     %[~,p] = ttest(kappa(:,1),kappa(:,2));
%     %p = ranksum(kappa(:,1),kappa(:,2));
%     text(r, 0.35,['p = ',num2str(p)]);
    
    
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
        
end

% figure(1)
% title('kappa')
% xticks([1:size(regions,1)]);
% xticklabels(xlabelstr)
% %axis([0 r+1 0 1]);
% figfile = [figDir,'PhaseLocking\',freq,'plPerf_kappa'];
% print('-dpdf', figfile);
% print('-djpeg', figfile);


figure
title('zrayl')
xlim([0 3])
xticks([0:size(regions,1)]);
xticklabels(xlabelstr)
%axis([0 r+1 -1 9]);
% if strcmp(freq,'beta')
%     ylim([0 3])
% end
% if strcmp(freq,'resp')
%     ylim([0 30])
% end
%figfile = [figDir,'PhaseLocking\',freq,'plPerf_zrayl'];
%print('-dpdf', figfile);
%print('-djpeg', figfile);
