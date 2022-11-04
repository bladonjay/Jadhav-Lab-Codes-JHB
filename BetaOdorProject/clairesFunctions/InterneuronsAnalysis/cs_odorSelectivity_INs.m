%% 
clear
[topDir, figDir] = cs_setPaths();

win = [0 1];
prewin=[-1 0];
binsize = 0.05;
% cs_odorSelectivity_v2(topDir, figDir, win,binsize)
saveout=0;


animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
winsize = win(2) + win(1);

%load([figDir, 'cmap_selectivity.mat']);
load('redToBlue');

%%

for r = 1:length(regions)
    region = regions{r};
    %find cells that meet selection criteria (pyramidal cells in area of
    %interest)
    selectivityAllCells = []; inds =[];
    i_SelectivityAllCells = [];
    allSI = [];
    
    psth_left_correct = [];
    psth_right_correct = [];
    psth_left_incorrect = [];
    psth_right_incorrect = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        %Take previously defined selective cells
        cellfilter = ['((strcmp($type, ''int'')) && (isequal($area,''',region, ...
            ''')) && (strcmp($selectivity, ''leftSelective'') || strcmp($selectivity, ''rightSelective'')))'];
        
        cells = evaluatefilter(cellinfo,cellfilter);
        
        noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        
        days = unique(cells(:,1));
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            %runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
            
            for c = 1:size(daycells,1)
                
                correctleftspikes = [];
                correctrightspikes = [];
                incorrectleftspikes = [];
                incorrectrightspikes = [];
                
                cell = daycells(c,:);
                runeps = cs_findGoodEpochs(cellinfo{day}, {'SI'},cell(2:3));
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                        epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                        
                        [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                        %triginds = sort([correct_left; correct_right]);
                        
                        trigs = odorTriggers{day}{epoch}.allTriggers;
                        %% Find Spikes on each trial
                        for t = 1:length(trigs)
                            trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
                            winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                            bins = (trigs(t)-win(1):binsize:trigs(t)+win(2));
                            binspikecount = histcounts(winspikes,bins);
                            
                            %classify into trial type
                            if ismember(t, correct_left)
                                correctleftspikes = [correctleftspikes; binspikecount];
                            elseif ismember(t, correct_right)
                                correctrightspikes = [correctrightspikes; binspikecount];
                            elseif ismember(t, incorrect_left)
                                incorrectleftspikes = [incorrectleftspikes; binspikecount];
                            elseif ismember(t, incorrect_right)
                                incorrectrightspikes = [incorrectrightspikes; binspikecount];
                            end
                            
                        end
                        
                    end
                    
                end
                
                cleftbinfr = mean(correctleftspikes,1)./binsize;
                crightbinfr = mean(correctrightspikes,1)./binsize;
                
                psth_left_correct = [psth_left_correct;cleftbinfr];
                psth_right_correct = [psth_right_correct;crightbinfr];
                
                seldat = (cleftbinfr - crightbinfr)./(cleftbinfr + crightbinfr);
                
                selectivityAllCells = [selectivityAllCells;seldat];
                
                inds = [inds; a, cell];
                
                
                ileftbinfr = mean(incorrectleftspikes,1)./binsize;
                if isempty(ileftbinfr)
                    ileftbinfr = zeros(1,length(bins)-1);
                end
                irightbinfr = mean(incorrectrightspikes,1)./binsize;
                if isempty(irightbinfr)
                    ibinfr = zeros(1,length(bins)-1);
                end
                
                
                if isempty(ileftbinfr) && isempty(irightbinfr)
                    i_SelectivityAllCells = [i_SelectivityAllCells; zeros(1,length(bins)-1)];
                else
                    
                    psth_left_incorrect = [psth_left_incorrect;ileftbinfr];
                    psth_right_incorrect = [psth_right_incorrect;irightbinfr];
                    
                    i_selectivity = (ileftbinfr - irightbinfr)./(ileftbinfr + irightbinfr);
                    i_SelectivityAllCells = [i_SelectivityAllCells; i_selectivity];
                end
                
                allSI = [allSI; cellinfo{cell(1)}{runeps(1)}{cell(2)}{cell(3)}.SI];
                
            end
            
        end
        
        
    end
    selectivityAllCells(isnan(selectivityAllCells))=0;
    
    bins = -win(1):binsize:win(2);
    trigbins = (bins >= 0.2 & bins < 1);
    mn = mean(selectivityAllCells(:,trigbins), 2);
    [corr,sortmean] = sort(mn,'descend');
    cellinds = inds(sortmean,:);
    
    
    %validinds = 1+((numpoints-1)/2):size(selectivityAllCells,2)-((numpoints-1)/2);
    %smoothed = zeros(length(validinds));
    newbins = bins(1)+binsize:binsize/5:bins(end);
    %1:0.25:length(selectivityAllCells(c,:));
    smoothed = zeros(size(selectivityAllCells,1),length(newbins));
    
    for c = 1:size(selectivityAllCells,1)
        
        %rather than smoothing, use interpolation instead.
        %smoothed(c,:) = smoothdata(selectivityAllCells(c,:),'gaussian',5);
        smoothed(c,:) = smoothdata(interp1(bins(2:end),selectivityAllCells(c,:),newbins),'gaussian',15);
        
    end
    
    
    smoothed = smoothed(sortmean,:);
    figure,
    set(gcf,'Position',[300 100 300 500]);
    imagesc(-(win(1))+binsize/5:binsize:win(2), length(selectivityAllCells):1,smoothed);
    colorbar
    %     [cmap]=buildcmap('rkg');
    colormap(redToBlue) %will use the output colormap
    colorbar('YTick', [-1 0 1]);
    caxis([-1 1])
    title([region, ' Correct Trials']);
    xlabel('Time from odor onset (s)');
    ylabel('Cell Number');
    
    %save
    if saveout==1
        figfile = [figDir,'Interneurons\OdorSelectivity_',region,'_Correct'];
        %saveas(gcf,figfile,'fig');
        print('-dpdf', figfile);
        print('-djpeg',figfile);
    end
    
    %Compare to Incorrect using same cell order
    
    
    
    %validinds = 1+((numpoints-1)/2):size(selectivityAllCells,2)-((numpoints-1)/2);
    
    
    i_SelectivityAllCells(isnan(i_SelectivityAllCells))=0;
    %     buffer = [repmat(i_SelectivityAllCells(:,1), [1, buffersize]), repmat(i_SelectivityAllCells(:,end), [1, buffersize])];
    %     i_SelectivityAllCells = [buffer(:,1:buffersize), i_SelectivityAllCells, buffer(:,buffersize+1:end)];
    %     i_smoothed = zeros(length(validinds));
    i_smoothed = zeros(size(i_SelectivityAllCells,1),length(newbins));
    
    
    for c = 1:size(i_SelectivityAllCells,1)
        %i_smoothed(c,:) = filter2(s,i_SelectivityAllCells(c,:),'valid');
        %i_smoothed(c,:) = smoothdata(i_SelectivityAllCells(c,:),'gaussian',5);
        i_smoothed(c,:) = smoothdata(interp1(bins(2:end),i_SelectivityAllCells(c,:),newbins),'gaussian',15);
    end
    i_smoothed = i_smoothed(sortmean,:);
    mn_i = mean(i_SelectivityAllCells(:,trigbins), 2);
    incorr = mn_i(sortmean);
    figure,
    set(gcf,'Position',[300 100 300 500]);
    imagesc([-(win(1)):binsize:win(2)], [length(i_SelectivityAllCells):1],i_smoothed);
    colorbar
    colorbar('YTick', [-1 0 1]);
    caxis([-1 1])
    %[cmap]=buildcmap('rkg');
    %colormap(cmap) %will use the output colormap
    colormap(redToBlue)
    title([region, ' Incorrect Trials']);
    xlabel('Time from odor onset (s)');
    ylabel('Cell Number');
    
    %save
    if saveout==1
        figfile = [figDir,'Interneurons\OdorSelectivity_',region,'_Incorrect'];
        %saveas(gcf,figfile,'fig');
        print('-dpdf', figfile);
        print('-djpeg',figfile);
    end
    
    %% Calculate correlation
    figure,
    plot(corr,incorr,'k.','MarkerSize',20)
    hold on
    fit = polyfit(corr, incorr,1);
    plot([-1 1], polyval(fit,[-1, 1]))
    xlabel('Correct Trial Selectiviy Index');
    ylabel('Incorrect Trial Selectivity Index');
    axis([-1 1 -1 1]);
    
    
    [CC,p] = corrcoef(corr, incorr);
    R = CC(1,2);
    p = p(1,2);
    
    txt = {['R = ',num2str(R)],['p = ' num2str(p)]};
    text(0.4,0.5,txt)
    if saveout==1
        figfile = [figDir,'Interneurons\OdorSelectivity_',region,'_SICorrelation'];
        %saveas(gcf,figfile,'fig');
        print('-dpdf', figfile);
        print('-djpeg',figfile);
    end
    
    %% save data
    
    %sort psth before saving
    psth_left_correct = psth_left_correct(sortmean,:);
    psth_right_correct = psth_right_correct(sortmean,:);
    psth_left_incorrect = psth_left_incorrect(sortmean,:);
    psth_right_incorrect = psth_right_incorrect(sortmean,:);
    
    
    selectivityData.psth_left_correct = psth_left_correct;
    selectivityData.psth_right_correct = psth_right_correct;
    selectivityData.psth_left_incorrect = psth_left_incorrect;
    selectivityData.psth_right_incorrect = psth_right_incorrect;
    selectivityData.SI_correct = smoothed;
    selectivityData.SI_incorrect = i_smoothed;
    selectivityData.SI_correct_mean = corr;
    selectivityData.SI_incorrect_mean = incorr;
    selectivityData.cellinds = cellinds;
    selectivityData.win = newbins;
    save([topDir,'AnalysesAcrossAnimals\selectivityData_IN_',region],'selectivityData');
    clear selectivityData
end


