%cs_phaseLocking_performance_v2
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

%Does this for all cells together, not split up into different eeg regions

clear
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};


%do all regions separately?
cellregions = {'CA1','PFC'};
eegregions = {'CA1','PFC'};

%get all region combinations
[A,B] = meshgrid(cellregions,eegregions);
c=cat(2,A',B');
regions = reshape(c,[],2);

for cr = 1:length(cellregions)
    cellregion = cellregions{cr};
    
    zrayl = [];
    kappa = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir =([topDir, animal,'Expt\',animal, '_direct\']);
        days = cs_getRunEpochs(animDir, animal, 'odorplace');
        days = unique(days(:,1));

        plfiles_c = dir([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-','*']);
        names_c = {plfiles_c.name};
        
        plfiles_i = dir([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_',cellregion,'-','*']);
        names_i = {plfiles_i.name};
        
        for d = days'
            daystr = getTwoDigitNumber(d);
            dayfiles_c = {plfiles_c(find(contains(names_c,daystr))).name};
            dayfiles_i = {plfiles_i(find(contains(names_i,daystr))).name};
            
            if isempty(dayfiles_c) || isempty(dayfiles_i)
                continue
            end
            
            allcells = [];
            for j = 1:length(dayfiles_c)
                dayfile_c = dayfiles_c{j};
%                 dayfile_i = dayfiles_i{j};
                if contains(dayfile_c, 'OB') %|| contains(dayfile_i, 'OB')
                    continue
                end
                
                load([animDir,'PhaseLocking\',dayfile_c]);
                eval(['phaselock_c = beta_phaselock',cellregion,';']);
                
%                 load([animDir,'PhaseLocking\Incorrect\',dayfile_i]);
%                 eval(['phaselock_i = beta_phaselock',cellregion,';']);
%                 
         
            %get phase locked cells. use cells that are phase locked on
            %correct trials. 
            cellfilt = '~isempty($sph) & $prayl <0.05';
            cells = evaluatefilter(phaselock_c,cellfilt);
            
            allcells = [allcells;cells];
            end
            
            
            
            if isempty(allcells)
                continue
            end
            cells = unique(allcells(:,[1 3 4]),'rows');
            
            %For each cell, get avg kappa and zrayl.
            for c = 1:size(cells,1)
                ind = cells(c,:);
                daystr = getTwoDigitNumber(ind(1));
                
                %% get avg
                k = [];
                z = [];
                for er = 1:length(eegregions)
                    eegregion = eegregions{er};
                    load([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-',eegregion,'_',daystr])
                    eval(['pl_c = beta_phaselock',cellregion,';']);
                    eps_c = cs_findGoodEpochs(pl_c{ind(1)},{'zrayl'},ind([2,3]));
                    
                    load([animDir,'PhaseLocking\Incorrect\',animal,'betaphaselock_',cellregion,'-',eegregion,'_incorrect_',daystr])
                    eval(['pl_i = beta_phaselock',cellregion,';']);
                    eps_i = cs_findGoodEpochs(pl_i{ind(1)},{'zrayl'},ind([2,3]));
                    
                    eps = intersect(eps_c, eps_i);
                    
                    
                    for ep = 1:length(eps)
                    epoch = eps(ep);
                    %get kappa and zrayl scores for comparison
                    k_ep = [pl_c{ind(1)}{epoch}{ind(2)}{ind(3)}.kappa, mean(pl_i{ind(1)}{epoch}{ind(2)}{ind(3)}.kappa)];
                    z_ep = [pl_c{ind(1)}{epoch}{ind(2)}{ind(3)}.zrayl, mean(pl_i{ind(1)}{epoch}{ind(2)}{ind(3)}.zrayl)];
                
                    k = [k;k_ep];
                    z = [z;z_ep];
 
                    end               
                end
                kappa = [kappa; mean(k,1)];
                zrayl = [zrayl; mean(z,1)];
            end
        end
    end
              
 
%     % Kappa figure
%     figure(f1)
%     k_std = std(kappa,1);
%     k_mean = mean(kappa,1);
%     errorbar(k_mean,k_std,'k.')
%     
%     figtitle = [cellregion,' kappa'];
%     k_std = std(kappa,1);
%     k_mean = mean(kappa,1);
%     errorbar(k_mean,k_std,'k.')
%     axis([0 3 0 (max(k_mean(:)) +2*(max(k_std(:))))])
%     xticks([1 2])
%     xticklabels({'Correct','Incorrect'})
%     [~,p] = ttest(kappa(:,1),kappa(:,2));
%     text(1.5, 3*(max(k_std(:))),['p = ',num2str(p)]);
%     title(figtitle)
%     figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
%     print('-dpdf', figfile);
%     print('-djpeg', figfile);
    
     % Zrayl figure
    figtitle = [cellregion,' zrayl'];
    z_std = std(zrayl,1);
    z_mean = mean(zrayl,1);
    errorbar(z_mean,z_std,'k.')
    axis([0 3 0 (max(z_mean(:)) +2*(max(z_std(:))))])
    xticks([1 2])
    xticklabels({'Correct','Incorrect'})
    [~,p] = ttest(zrayl(:,1),zrayl(:,2));
    text(1.5, 3*(max(z_std(:))),['p = ',num2str(p)]);
    title(figtitle)
    figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
end