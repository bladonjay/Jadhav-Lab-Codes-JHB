%cs_phaseLocking_performance
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

%uses all cells that are phase locked (not just to one region)
%but only takes rayleigh Z from local region

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
eegregions = {'CA1','PFC'};
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
    
    
    %get phase locked cells. use cells that are phase locked on
        %correct trials.
        allplcells = [];
        for re = 1:length(eegregions)
            eegr = eegregions{re};
            load([topDir, 'AnalysesAcrossAnimals\PhaseLocking\plCells_',freq,'_',cellregion,'-',eegr,'.mat']);
            allplcells = [allplcells; plcells];
        end
        cells = unique(allplcells,'rows');
                
    for a = 1:length(animals)
        animal = animals{a};
        animDir =([topDir, animal,'Expt\',animal, '_direct\']);
        
        animCells = (cells(cells(:,1) == a,:));
        if isempty(animCells)
            continue
        end
        
        load([animDir,'PhaseLocking\',animal,'phaselock_',freq,'_',cellregion,'-',eegregion]);
        phaselock_correct = phaselock;
        try
            load([animDir,'PhaseLocking\Incorrect\',animal,'phaselock_',freq,'_incorrect_',cellregion,'-',eegregion]);
            phaselock_incorrect = phaselock;
        catch
            continue
        end

        for c = 1:size(animCells,1)
            ind = animCells(c,[2:4]);
            allcells = [allcells;a,ind];
            %
%             k = phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.kappa;
%             kdist = phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.kappa_dist;
            
            z_i = phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.zrayl;
            z_c = mean(phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.zrayl_dist);
            
            k_i = phaselock_incorrect{ind(1)}{1}{ind(2)}{ind(3)}.kappa;
            k_c = mean(phaselock_correct{ind(1)}{1}{ind(2)}{ind(3)}.kappa);
            
            
            
%             kappa = [kappa;k];
%             kappa_dist = [kappa_dist, kdist];
            
            zrayl = [zrayl;z_c,z_i];
            kappa = [kappa;k_c,k_i];
            %zrayl_dist = [zrayl_dist, zdist];
            
            
            
            
%             if isempty(kappa) || isempty(zrayl)
%                 keyboard
%             end
        end
        
    end
   
%     
    
    % Zrayl figure
    figure(1)
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
    figure(2)
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
        
end

% figure(1)
% title('kappa')
% xticks([1:size(regions,1)]);
% xticklabels(xlabelstr)
% %axis([0 r+1 0 1]);
% figfile = [figDir,'PhaseLocking\',freq,'plPerf_kappa'];
% print('-dpdf', figfile);
% print('-djpeg', figfile);


figure(1)
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
figfile = [figDir,'PhaseLocking\',freq,'plPerf_zrayl'];
print('-dpdf', figfile);
print('-djpeg', figfile);
