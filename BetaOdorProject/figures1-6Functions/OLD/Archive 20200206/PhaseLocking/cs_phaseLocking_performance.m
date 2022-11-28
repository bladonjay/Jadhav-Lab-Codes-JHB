%cs_phaseLocking_performance
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

clear
close all
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

freq = 'resp';

%do all regions separately?
cellregions = {'CA1','PFC'};
eegregions = {'CA1','PFC','OB'};
%get all region combinations
[A,B] = meshgrid(cellregions,eegregions);
c=cat(2,A',B');
regions = reshape(c,[],2);
xlabelstr = [];
for r = 1:size(regions,1)
    cellregion = regions{r,1};
    eegregion = regions{r,2};
    
    zrayl = [];
    kappa = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir =([topDir, animal,'Expt\',animal, '_direct\']);
        plfiles = dir([animDir,'PhaseLocking\',animal,freq,'phaselock_',cellregion,'-',eegregion,'_*']);
        for f = 1:length(plfiles)
            load([animDir,'PhaseLocking\',plfiles(f).name])
            eval(['phaselock_correct = ',freq,'_phaselock',cellregion,';']);
            
            daystr = plfiles(f).name(end-5:end-4);
            try
                load([animDir,'PhaseLocking\Incorrect\',animal,freq,'phaselock_',cellregion,'-',eegregion,'_incorrect_',daystr]);
                eval(['phaselock_incorrect = ',freq,'_phaselock',cellregion,';']);
            catch
                continue
            end
            
            %get phase locked cells. use cells that are phase locked on
            %correct trials.
            cellfilt = '~isempty($sph) & $prayl <0.05';
            cells = evaluatefilter(phaselock_correct,cellfilt);
            
            if isempty(cells)
                continue
            end
            
            cells = unique(cells(:,[1 3 4]),'rows');
            
            
            
            for c = 1:size(cells,1)
                ind = cells(c,:);
                eps_c = cs_findGoodEpochs(phaselock_correct{ind(1)},{'sph'},ind([2,3]));
                eps_i = cs_findGoodEpochs(phaselock_incorrect{ind(1)},{'sph'},ind([2,3]));
                eps = intersect(eps_c, eps_i); %make sure cell was present in both correct and incorrect trials
                if isempty(eps)
                    continue
                end
                
%                 if isempty(eps_i) && ~isempty(eps_c)
%                     kappa = [kappa; phaselock_correct{ind(1)}{eps_c}{ind(2)}{ind(3)}.kappa, 0];
%                     zrayl = [zrayl; phaselock_correct{ind(1)}{eps_c}{ind(2)}{ind(3)}.zrayl, 0];
%                 end
                
                %May need to change how high beta periods are defined. Beta
                %power power will be lower on incorrect trials to begin
                %with...
                k = []; z = [];
                for ep = 1:length(eps)
                    epoch = eps(ep);
                    %get kappa and zrayl scores for comparison
                    k_ep = [phaselock_correct{ind(1)}{epoch}{ind(2)}{ind(3)}.kappa, mean(phaselock_incorrect{ind(1)}{epoch}{ind(2)}{ind(3)}.kappa)];
                    z_ep = [phaselock_correct{ind(1)}{epoch}{ind(2)}{ind(3)}.zrayl, mean(phaselock_incorrect{ind(1)}{epoch}{ind(2)}{ind(3)}.zrayl)];
                    
                    k = [k;k_ep];
                    z = [z;z_ep];
                    
                    if isempty(k) || isempty(z)
                        keyboard
                    end
                end
                kappa = [kappa;mean(k,1)];
                zrayl = [zrayl;mean(z,1)];
                
                if isempty(kappa) || isempty(zrayl)
                    keyboard
                end
            end
        end
        
    end
    
    if isempty(kappa) || isempty(zrayl)
        continue
    end
    
    % Kappa figure
    figure(1)
    hold on
    k_std = std(kappa,1);
    k_mean = mean(kappa,1);
    errorbar([r, r+0.25], k_mean,k_std,'k.')
    xlabelstr = [xlabelstr,{[cellregion,'-',eegregion]}];
    %     figtitle = [cellregion,'-',eegregion,' kappa'];
    %     k_std = std(kappa,1);
    %     k_mean = mean(kappa,1);
    %     errorbar(k_mean,k_std,'k.')
    %     axis([0 3 0 (max(k_mean(:)) +2*(max(k_std(:))))])
    %     xticks([1 2])
    %xticklabels({'Correct','Incorrect'})
    [~,p] = ttest(kappa(:,1),kappa(:,2));
    text(r, 1,['p = ',num2str(p)]);
    
    
    % Zrayl figure
    figure(2)
    hold on
    %     figtitle = [cellregion,'-',eegregion,' zrayl'];
    z_std = stderr(zrayl);
    z_mean = mean(zrayl,1);
    errorbar([r, r+0.25], z_mean,z_std,'k.')
    %     axis([0 3 0 (max(z_mean(:)) +2*(max(z_std(:))))])
    %     xticks([1 2])
    %     xticklabels({'Correct','Incorrect'})
    [~,p] = ttest(zrayl(:,1),zrayl(:,2));
    text(r,1,['p = ',num2str(p)]);
    %     title(figtitle)
    %     figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
    %     print('-dpdf', figfile);
    %     print('-djpeg', figfile);
end

figure(1)
title('kappa')
xticks([1:size(regions,1)]);
xticklabels(xlabelstr)
%axis([0 r+1 0 1]);
figfile = [figDir,'PhaseLocking\',freq,'plPerf_kappa'];
print('-dpdf', figfile);
print('-djpeg', figfile);


figure(2)
title('zrayl')
xticks([1:size(regions,1)]);
xticklabels(xlabelstr)
%axis([0 r+1 -1 9]);
xlim([0 5])
figfile = [figDir,'PhaseLocking\',freq,'plPerf_zrayl'];
print('-dpdf', figfile);
print('-djpeg', figfile);
