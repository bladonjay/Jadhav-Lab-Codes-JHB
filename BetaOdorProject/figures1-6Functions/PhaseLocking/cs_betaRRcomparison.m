%cs_phaseLocking_performance
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

clear
close all
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};


regions = {'CA1','CA1';'PFC','PFC'};
xlabelstr = [];
for r = 1:size(regions,1)
    cellregion = regions{r,1};
    eegregion = regions{r,2};
    
    zrayl = [];
    kappa = [];
    R_all = [];
    r_dist = [];
    allcells = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir =([topDir, animal,'Expt\',animal, '_direct\']);
        
        load([animDir,'PhaseLocking\',animal,'phaselock_beta_',cellregion,'-',eegregion]);
        beta_phaselock = phaselock;
        load([animDir,'PhaseLocking\',animal,'phaselock_resp_',cellregion,'-',eegregion]);
        resp_phaselock = phaselock;
        
        cellfilt = '~isempty($sph)';
        
        cells = evaluatefilter(beta_phaselock,cellfilt);
        
        if isempty(cells)
            continue
        end
        
        cells = unique(cells(:,[1 3 4]),'rows');
        
        
        
        for c = 1:size(cells,1)
            ind = cells(c,:);
            allcells = [allcells;a,ind];
            
            k_b = beta_phaselock{ind(1)}{1}{ind(2)}{ind(3)}.kappa;
            k_r = resp_phaselock{ind(1)}{1}{ind(2)}{ind(3)}.kappa;
            kappa = [kappa;[k_b,k_r]];
            
            z_b = beta_phaselock{ind(1)}{1}{ind(2)}{ind(3)}.zrayl;
            z_r = resp_phaselock{ind(1)}{1}{ind(2)}{ind(3)}.zrayl;
            zrayl = [zrayl;[z_b,z_r]];
            
        end
        
        %
        
    end
    
    
    figure, hold on
    p = ranksum(kappa(:,1),kappa(:,2));
    sem = stderr(kappa);
    mn = mean(kappa);
    bar([1 2],mn)
    errorbar([1,2],mn,sem);
    text(1,mn(1),['p = ',num2str(p)]);
    xticks([1 2])
    xticklabels({'beta','RR'});
    ylabel('kappa')
    xlim([0 3])
    ylim([0 inf])
    figfile = [figDir,'PhaseLocking\betaRRcomp_kappa_',cellregion];
    %print('-dpdf', figfile);
    %print('-djpeg', figfile);
    
    figure, hold on
    p = ranksum(zrayl(:,1),zrayl(:,2));
    sem = stderr(zrayl);
    mn = mean(zrayl);
    bar([1 2],mn)
    errorbar([1,2],mn,sem);
    text(1,mn(1),['p = ',num2str(p)]);
      xticks([1 2])
    xticklabels({'beta','RR'});
    ylabel('Rayleigh Z')
    xlim([0 3])
    ylim([0 inf])
    figfile = [figDir,'PhaseLocking\betaRRcomp_zrayl_',cellregion];
    %print('-dpdf', figfile);
    %print('-djpeg', figfile);
    
end

