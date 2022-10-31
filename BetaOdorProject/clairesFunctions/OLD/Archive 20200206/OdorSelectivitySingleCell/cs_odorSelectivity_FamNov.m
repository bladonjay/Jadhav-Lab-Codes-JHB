%% Params
clear
[topDir, figDir] = cs_setPaths();

win = [0.5 1];
binsize = 0.05;
% cs_odorSelectivity_v2(topDir, figDir, win,binsize)



animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
winsize = win(2) + win(1);

load([figDir, 'cmap_selectivity.mat']);

%%

for r = 1:length(regions)
    region = regions{r};
    %find cells that meet selection criteria (pyramidal cells in area of
    %interest)
    selectivity_familiar = []; inds =[];
    selectivity_novel = [];
    allSI = [];
    load([topDir,'AnalysesAcrossAnimals\npCells_',region]);
    for a = 5%1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        cells = npCells(npCells(:,1) == a,2:4);
        
        load([animDir,animal,'cellinfo.mat'])
        
        %         %Take previously defined selective cells
        %         cellfilter = ['((strcmp($type, ''pyr'')) && (isequal($area,''',region, ...
        %             ''')) && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
        %
        %         cells = evaluatefilter(cellinfo,cellfilter);
        
        %noeps = cells(:,[1 3 4]);
        %cells = unique(noeps,'rows');
        
        days = unique(cells(:,1));
        for d = 2:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            spikes = loaddatastruct(animDir,animal,'spikes',day);
            odorTriggers = loaddatastruct(animDir,animal,'odorTriggers',day);
            %runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
            familiarEps = cs_getRunEpochs(animDir,animal,'odorplace',day);
            novel1Eps = cs_getRunEpochs(animDir,animal,'novelodor',day);
            novel2Eps = cs_getRunEpochs(animDir,animal,'novelodor2',day);
            
            novelEps = [novel1Eps;novel2Eps];
            
            for c = 1:size(daycells,1)
                
                leftspikes = [];
                rightspikes = [];
                %                 incorrectleftspikes = [];
                %                 incorrectrightspikes = [];
                
                cell = daycells(c,:);
                %runeps = cs_findGoodEpochs(cellinfo{day}, 'SI',cell(2:3));
                
                for ep = 1:size(familiarEps,1)
                    epoch = familiarEps(ep,2);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                        epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                        
                        [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                        %triginds = sort([correct_left; correct_right]);
                        
                        trigs = odorTriggers{day}{epoch}.allTriggers(sort([correct_left;correct_right]));
                        %% Find Spikes on each trial
                        for t = 1:length(trigs)
                            trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
                            winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                            bins = (trigs(t)-win(1):binsize:trigs(t)+win(2));
                            binspikecount = histcounts(winspikes,bins);
                            
                            %classify into trial type
                            if ismember(t, correct_left)
                                leftspikes = [leftspikes; binspikecount];
                            elseif ismember(t, correct_right)
                                rightspikes = [rightspikes; binspikecount];
                                %                             elseif ismember(t, incorrect_left)
                                %                                 incorrectleftspikes = [incorrectleftspikes; binspikecount];
                                %                             elseif ismember(t, incorrect_right)
                                %                                 incorrectrightspikes = [incorrectrightspikes; binspikecount];
                            end
                            
                        end
                        
                    end
                    
                end
                
                leftbinfr = mean(leftspikes,1)./binsize;
                rightbinfr = mean(rightspikes,1)./binsize;
                
                selectivity = (leftbinfr - rightbinfr)./(leftbinfr + rightbinfr);
                
                selectivity_familiar = [selectivity_familiar;selectivity];

                %% Calculate SI for novel odor trials
                leftspikes = [];
                rightspikes = [];
                
                for ep = 1:size(novelEps,1)
                    epoch = novelEps(ep,2);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                        epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                        
                        [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                        %triginds = sort([correct_left; correct_right]);
                        
                        trigs = odorTriggers{day}{epoch}.allTriggers(sort([correct_left;correct_right]));
                        %% Find Spikes on each trial
                        for t = 1:length(trigs)
                            trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
                            winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                            bins = (trigs(t)-win(1):binsize:trigs(t)+win(2));
                            binspikecount = histcounts(winspikes,bins);
                            
                            %classify into trial type
                            if ismember(t, correct_left)
                                leftspikes = [leftspikes; binspikecount];
                            elseif ismember(t, correct_right)
                                rightspikes = [rightspikes; binspikecount];
                                %                             elseif ismember(t, incorrect_left)
                                %                                 incorrectleftspikes = [incorrectleftspikes; binspikecount];
                                %                             elseif ismember(t, incorrect_right)
                                %                                 incorrectrightspikes = [incorrectrightspikes; binspikecount];
                            end
                            
                        end
                        
                    end
                    
                end
                leftbinfr = mean(leftspikes,1)./binsize;
                rightbinfr = mean(rightspikes,1)./binsize;
                
                selectivity = (leftbinfr - rightbinfr)./(leftbinfr + rightbinfr);
                if isempty(selectivity)
                    selectivity = zeros(1,length(bins)-1);
                end
                
                selectivity_novel = [selectivity_novel;selectivity];
                
                %inds = [inds; a, cell];
                %allSI = [allSI; cellinfo{cell(1)}{runeps(1)}{cell(2)}{cell(3)}.SI];
                
            end
            
        end
        
        
    end
    selectivity_familiar(isnan(selectivity_familiar))=0;
    
    bins = -win(1):binsize:win(2);
    trigbins = (bins >= 0.2 & bins < 1);
    mn = mean(selectivity_familiar(:,trigbins), 2);
    [fam,sortmean] = sort(mn,'descend');
    
    smoothed = zeros(size(selectivity_familiar,1),size(selectivity_familiar,2));
    
    for c = 1:size(selectivity_familiar,1)
        smoothed(c,:) = smoothdata(selectivity_familiar(c,:),'gaussian',5);
    end
    
    smoothed = smoothed(sortmean,:);
    figure,
    set(gcf,'Position',[300 100 300 500]);
    imagesc([-(win(1)):binsize:win(2)], [length(selectivity_familiar):1],smoothed);
    colorbar
    %     [cmap]=buildcmap('rkg');
    colormap(cmap) %will use the output colormap
    colorbar('YTick', [-1 0 1]);
    caxis([-1 1])
    title([region, ' Familiar Odors']);
    xlabel('Time from odor onset (s)');
    ylabel('Cell Number');
    
    %save
    figfile = [figDir,'IgarashiReplication\4a_OdorSelectivity_',region,'_Familiar'];
    saveas(gcf,figfile,'fig');
    print('-dpdf', figfile);
    print('-djpeg',figfile);
    
    
    %Compare to Novel using same cell order
   
    
    selectivity_novel(isnan(selectivity_novel))=0;
   
    i_smoothed = zeros(size(selectivity_novel,1),size(selectivity_novel,2));
    
    
    for c = 1:size(selectivity_novel,1)
        i_smoothed(c,:) = smoothdata(selectivity_novel(c,:),'gaussian',5);
    end
    i_smoothed = i_smoothed(sortmean,:);
    mn_i = mean(selectivity_novel(:,trigbins), 2);
    nov = mn_i(sortmean);
    figure,
    set(gcf,'Position',[300 100 300 500]);
    imagesc([-(win(1)):binsize:win(2)], [length(selectivity_novel):1],i_smoothed);
    colorbar
    colorbar('YTick', [-1 0 1]);
    caxis([-1 1])
    %[cmap]=buildcmap('rkg');
    %colormap(cmap) %will use the output colormap
    colormap(cmap)
    title([region, ' Novel Odors']);
    xlabel('Time from odor onset (s)');
    ylabel('Cell Number');
    
    %save
    figfile = [figDir,'IgarashiReplication\4a_OdorSelectivity_',region,'_Novel'];
    %saveas(gcf,figfile,'fig');
    print('-dpdf', figfile);
    print('-djpeg',figfile);
    
    
    %% Calculate correlation
    figure,
    plot(fam,nov,'k.','MarkerSize',20)
    hold on
    fit = polyfit(fam, nov,1);
    plot([-1 1], polyval(fit,[-1, 1]))
    xlabel('Familiar Odor Selectiviy Index');
    ylabel('Novel Odor Selectivity Index');
    axis([-1 1 -1 1]);
    
    
    [CC,p] = corrcoef(fam, nov);
    R = CC(1,2)
    p = p(1,2)
    
    txt = {['R = ',num2str(R)],['p = ' num2str(p)]};
    text(0.4,0.5,txt)
    
    %save
    figfile = [figDir,'IgarashiReplication\4a_OdorSelectivity_',region,'_FamNovCorrelation'];
    %saveas(gcf,figfile,'fig');
    print('-dpdf', figfile);
    print('-djpeg',figfile);
end


