%% Params
[topDir, figDir] = cs_setPaths();

win = [0.5 1];
binsize = 0.05;
% cs_odorSelectivity_v2(topDir, figDir, win,binsize)



animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42'};
regions = {'CA1','PFC'};
winsize = win(2) + win(1);

load([figDir, 'cmap_selectivity.mat']);

%%
for r = 1:length(regions)
region = regions{r};
    %find cells that meet selection criteria (pyramidal cells in area of
    %interest)
    selectivityAllCells = []; inds =[];
    i_SelectivityAllCells = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        %Take previously defined selective cells
        cellfilter = ['((strcmp($type, ''pyr'')) && (isequal($area,''',region, ... 
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
            runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
           
            for c = 1:size(daycells,1)

                correctleftspikes = []; 
                correctrightspikes = [];
                incorrectleftspikes = [];
                incorrectrightspikes = [];
                
                cell = daycells(c,:);
                
                for ep = 1:length(runeps)
                    epoch = runeps(ep);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
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
                    
                    selectivity = (cleftbinfr - crightbinfr)./(cleftbinfr + crightbinfr);
                    
                    selectivityAllCells = [selectivityAllCells;selectivity];
                    
                    inds = [inds; a, cell];
                    
                    if isempty(incorrectleftspikes)
                        incorrectleftspikes = 0;
                    end
                    if isempty(incorrectrightspikes)
                        incorrectrightspikes = 0;
                    end
                    
                    ileftbinfr = mean(incorrectleftspikes,1)./binsize;
                    irightbinfr = mean(incorrectrightspikes,1)./binsize;
                    
                    i_selectivity = (ileftbinfr - irightbinfr)./(ileftbinfr + irightbinfr);
                    i_SelectivityAllCells = [i_SelectivityAllCells; i_selectivity]; 
                    
                %end
                
            end
            
        end
       

    end
    selectivityAllCells(isnan(selectivityAllCells))=0;
    
    bins = -win(1):binsize:win(2);
    trigbins = (bins >= 0.2 & bins < 0.8); 
    mn = mean(selectivityAllCells(:,trigbins), 2);
    [~,sortmean] = sort(mn,'descend');
    
    %std = 2; numpoints = 5; %number of points making up gaussian
    %s = gaussian(std,numpoints);
   % buffersize = (numpoints-1)/2;
    %buffer = [repmat(selectivityAllCells(:,1), [1, buffersize]), repmat(selectivityAllCells(:,end), [1, buffersize])];
    %selectivityAllCells = [buffer(:,1:buffersize), selectivityAllCells, buffer(:,buffersize+1:end)];
    
    %validinds = 1+((numpoints-1)/2):size(selectivityAllCells,2)-((numpoints-1)/2);
    %smoothed = zeros(length(validinds));
    smoothed = zeros(size(selectivityAllCells,1),size(selectivityAllCells,2));
    
    for c = 1:size(selectivityAllCells,1)
         %smoothed(c,:) = filter2(s,selectivityAllCells(c,:),'valid'); 
         smoothed(c,:) = smoothdata(selectivityAllCells(c,:),'gaussian',5); 
    end
    smoothed = smoothed(sortmean,:);
    figure,
    set(gcf,'Position',[300 100 300 500]);
    imagesc([-(win(1)):binsize:win(2)], [length(selectivityAllCells):1],smoothed);
    colorbar
%     [cmap]=buildcmap('rkg'); 
    colormap(cmap) %will use the output colormap
    colorbar('YTick', [-1 0 1]);
    caxis([-1 1])
    title([region, ' Correct Trials']);
    xlabel('Time from odor onset (s)');
    ylabel('Cell Number');
    
    %save
       figfile = [figDir,'IgarashiReplication\4a_OdorSelectivity_',region,'_Correct'];
       saveas(gcf,figfile,'fig');
       print('-dpdf', figfile);
       print('-djpeg',figfile);
    
    
    %Compare to Incorrect using same cell order
    
    
    
    %validinds = 1+((numpoints-1)/2):size(selectivityAllCells,2)-((numpoints-1)/2);
    
    
    i_SelectivityAllCells(isnan(i_SelectivityAllCells))=0;
%     buffer = [repmat(i_SelectivityAllCells(:,1), [1, buffersize]), repmat(i_SelectivityAllCells(:,end), [1, buffersize])];
%     i_SelectivityAllCells = [buffer(:,1:buffersize), i_SelectivityAllCells, buffer(:,buffersize+1:end)]; 
%     i_smoothed = zeros(length(validinds));
    i_smoothed = zeros(size(i_SelectivityAllCells,1),size(i_SelectivityAllCells,2));
    
    
    for c = 1:size(i_SelectivityAllCells,1)
         %i_smoothed(c,:) = filter2(s,i_SelectivityAllCells(c,:),'valid'); 
         i_smoothed(c,:) = smoothdata(i_SelectivityAllCells(c,:),'gaussian',5); 
    end
    i_smoothed = i_smoothed(sortmean,:);
    figure,
    set(gcf,'Position',[300 100 300 500]);
    imagesc([-(win(1)):binsize:win(2)], [length(i_SelectivityAllCells):1],i_smoothed);
    colorbar
    colorbar('YTick', [-1 0 1]);
    caxis([-1 1])
    %[cmap]=buildcmap('rkg'); 
    %colormap(cmap) %will use the output colormap
    colormap(cmap)
    title([region, ' Incorrect Trials']);
    xlabel('Time from odor onset (s)');
    ylabel('Cell Number');
    
    %save
       figfile = [figDir,'IgarashiReplication\4a_OdorSelectivity_',region,'_Incorrect'];
       %saveas(gcf,figfile,'fig');
       print('-dpdf', figfile);
       print('-djpeg',figfile);
    
end


