%cs_phaseLocking_performance_v2
%compare rayleigh z and kappa values for correct vs incorrect trials
%already have phase locking for all cells on correct/incorrect trials, just
%load files and get distributions

%Does this for all cells together, not split up into different eeg regions

clear
[topDir, figDir] = cs_setPaths();
animals = {'CS39','CS41','CS42','CS44'};


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
        dm1 = cs_getRunEpochs(animDir, animal, 'novelodor');
        dm2 = cs_getRunEpochs(animDir, animal, 'novelodor2');
        daymatrix = sortrows([dm1;dm2],1);
        
        days = unique(daymatrix(:,1));

        plfiles_pre = dir([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-','*']);
        names_pre = {plfiles_pre.name};
        
        plfiles_post = dir([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-','*']);
        names_post = {plfiles_post.name};
        
        for d = days'
            daystr = getTwoDigitNumber(d);
            dayfiles_pre = {plfiles_pre(find(contains(names_pre,['prelearn_',daystr]))).name};
            dayfiles_post = {plfiles_post(find(contains(names_post,['postlearn_',daystr]))).name};
            
            if isempty(dayfiles_pre) || isempty(dayfiles_post)
                continue
            end
            
            precells = [];
            for j = 1:length(dayfiles_pre)
                dayfile_pre = dayfiles_pre{j};
%                 dayfile_i = dayfiles_i{j};
                if contains(dayfile_pre, 'OB') %|| contains(dayfile_i, 'OB')
                    continue
                end
                
                load([animDir,'PhaseLocking\',dayfile_pre]);
                eval(['phaselock_pre = beta_phaselock',cellregion,';']);
                cellfilt = '~isempty($sph) & $prayl <0.05';
                cells = evaluatefilter(phaselock_pre,cellfilt);
            
                precells = [precells;cells];
            end
                
            cells = unique(precells(:,[1 3 4]),'rows');
            
            if isempty(cells)
                continue
            end

            %For each cell, get avg kappa and zrayl.
            for c = 1:size(cells,1)
                ind = cells(c,:);
                daystr = getTwoDigitNumber(ind(1));
                
                %% get avg
                k = [];
                z = [];
                for er = 1:length(eegregions)
                    eegregion = eegregions{er};
                    load([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-',eegregion,'_prelearn_',daystr])
                    eval(['pl_pre = beta_phaselock',cellregion,';']);
                    eps = cs_findGoodEpochs(pl_pre{ind(1)},{'sph'},ind([2,3]));
                    
%                     load([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-',eegregion,'_postlearn_',daystr])
%                     eval(['pl_post = beta_phaselock',cellregion,';']);
%                     eps_post = cs_findGoodEpochs(pl_post{ind(1)},{'sph'},ind([2,3]));
                    
%                     eps = intersect(eps_pre, eps_post);
                    
                    
                    for ep = 1:length(eps)
                    epoch = eps(ep);
                    %get kappa and zrayl scores for comparison
                    k_ep = pl_pre{ind(1)}{epoch}{ind(2)}{ind(3)}.kappa;
                    z_ep = pl_pre{ind(1)}{epoch}{ind(2)}{ind(3)}.zrayl;
                
                    k = [k;k_ep];
                    z = [z;z_ep];
 
                    end               
                end
                kappa = [kappa; mean(k,1)];
                zrayl = [zrayl; mean(z,1)];
            end
            
            
        end
    end
    k_pre = kappa;
    z_pre = zrayl;
    
    kappa = [];
    zrayl = [];
    %% POST learning
    for a = 1:length(animals)
        animal = animals{a};
        animDir =([topDir, animal,'Expt\',animal, '_direct\']);
        dm1 = cs_getRunEpochs(animDir, animal, 'novelodor');
        dm2 = cs_getRunEpochs(animDir, animal, 'novelodor2');
        daymatrix = sortrows([dm1;dm2],1);
        
        days = unique(daymatrix(:,1));
        
        plfiles_post = dir([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-','*']);
        names_post = {plfiles_post.name};
        
        for d = days'
            daystr = getTwoDigitNumber(d);
            dayfiles_post = {plfiles_post(find(contains(names_post,['postlearn_',daystr]))).name};
            
            if isempty(dayfiles_post)
                continue
            end
            
            postcells = [];
            for j = 1:length(dayfiles_post)
                dayfile_post = dayfiles_post{j};
%                 dayfile_i = dayfiles_i{j};
                if contains(dayfile_post, 'OB') %|| contains(dayfile_i, 'OB')
                    continue
                end
                
                load([animDir,'PhaseLocking\',dayfile_post]);
                eval(['phaselock_post = beta_phaselock',cellregion,';']);
                cellfilt = '~isempty($sph) & $prayl <0.05';
                cells = evaluatefilter(phaselock_post,cellfilt);
            
                postcells = [postcells;cells];
            end
                
            cells = unique(postcells(:,[1 3 4]),'rows');
            
            if isempty(cells)
                continue
            end

            %For each cell, get avg kappa and zrayl.
            for c = 1:size(cells,1)
                ind = cells(c,:);
                daystr = getTwoDigitNumber(ind(1));
                
                %% get avg
                k = [];
                z = [];
                for er = 1:length(eegregions)
                    eegregion = eegregions{er};
%                     load([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-',eegregion,'_prelearn_',daystr])
%                     eval(['pl_pre = beta_phaselock',cellregion,';']);
%                     eps_pre = cs_findGoodEpochs(pl_pre{ind(1)},{'sph'},ind([2,3]));
                    
                    load([animDir,'PhaseLocking\',animal,'betaphaselock_',cellregion,'-',eegregion,'_postlearn_',daystr])
                    eval(['pl_post = beta_phaselock',cellregion,';']);
                    eps = cs_findGoodEpochs(pl_post{ind(1)},{'sph'},ind([2,3]));

                    for ep = 1:length(eps)
                    epoch = eps(ep);
                    %get kappa and zrayl scores for comparison
                    k_ep = [pl_post{ind(1)}{epoch}{ind(2)}{ind(3)}.kappa];
                    z_ep = [pl_post{ind(1)}{epoch}{ind(2)}{ind(3)}.zrayl];
                
                    k = [k;k_ep];
                    z = [z;z_ep];
 
                    end               
                end
                kappa = [kappa; mean(k,1)];
                zrayl = [zrayl; mean(z,1)];
            end
        end
    end
    k_post = kappa;
    z_post = zrayl;
              
 
%     % Kappa figure
%     figtitle = [cellregion,' learning kappa'];
%     k_std = std(kappa,1);
%     k_mean = mean(kappa,1);
%     errorbar(k_mean,k_std,'k.')
%     axis([0 3 0 (max(k_mean(:)) +2*(max(k_std(:))))])
%     xticks([1 2])
%     xticklabels({'Prelearn','Postlearn'})
%     [~,p] = ttest(kappa(:,1),kappa(:,2));
%     text(1.5, 3*(max(k_std(:))),['p = ',num2str(p)]);
%     title(figtitle)
%     figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
%     print('-dpdf', figfile);
%     print('-djpeg', figfile);
    
     % Zrayl figure
    figtitle = [cellregion,' learning zrayl'];
    z_std = [std(z_pre) std(z_post)];
    z_mean = [mean(z_pre) mean(z_post)];
    errorbar(z_mean,z_std,'k.')
    axis([0 3 0 (max(z_mean(:)) +2*(max(z_std(:))))])
    xticks([1 2])
    xticklabels({'Prelearn','Postlearn'})
    [~,p] = ttest2(z_post,z_pre);
    text(1.5, 3*(max(z_std(:))),['p = ',num2str(p)]);
    title(figtitle)
    figfile = [figDir,'PhaseLocking\plPerf ',figtitle];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
end
