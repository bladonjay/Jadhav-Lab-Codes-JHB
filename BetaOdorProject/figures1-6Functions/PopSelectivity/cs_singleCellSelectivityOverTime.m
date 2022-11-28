%For each selective cell, define an average timepoint at which PDI is
%significant.
function cs_singleCellPDI(win, binsize)

% win = [.2 1];
% binsize = 0.05;

%% Params
animals = {'CS31','CS33','CS34','CS35'};

[topDir] = cs_setPaths();
regions = {'CA1','PFC'};

iterations = 1000;

plot = 0;
%%
for r = 1:length(regions)
    region = regions{r};
    
    %     load selective cells only
    load([topDir, 'AnalysesAcrossAnimals\selectiveCells_', region, '.mat'])
    cells = selectivecells;
    
    alltrigspikes = {};
    
    prevanimalday = [];
    for j = 1:size(cells,1)
        animaldayind = cells(j,[1 2]);
        
        
        if ~isequal(animaldayind, prevanimalday)
            animal = animals{cells(j,1)};
            day = cells(j,2);
            
            animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
            
            daystr = getTwoDigitNumber(day);
            load([animDir,animal,'spikes',daystr,'.mat'])
            load([animDir,animal,'odorTriggers',daystr,'.mat'])
            runeps = find(~cellfun(@isempty,odorTriggers{day}));
            
        end
        
        prevanimalday = animaldayind;
        
        cell = cells(j,[2 3 4]);
        
        epochleft = []; epochright = [];
        %get trials across day, matrix of ones and zeros. Then create
        %matrix with shuffled indicies, so its the same each time for each
        %cell on each day.
        totalleft = 0; totalright = 0;
        for ep = 1:length(runeps)
            epoch = runeps(ep);
            
            [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
            totalleft = totalleft + length(correct_left);
            totalright = totalright + length(correct_right);
            
            trigs = odorTriggers{day}{epoch}.allTriggers;
            
            lefttrigs = trigs(correct_left);
            righttrigs = trigs(correct_right);
            
            leftspikes = [];
            rightspikes = [];
            if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                
                
                for t = 1:length(lefttrigs)
                    trigwin = [lefttrigs(t)-win(1), lefttrigs(t)+win(2)];
                    winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                    bins = (lefttrigs(t)-win(1):binsize:lefttrigs(t)+win(2));
                    binspikes = histcounts(winspikes',bins);
                    leftspikes = [leftspikes; binspikes];
                    
                end
                
                for t = 1:length(righttrigs)
                    trigwin = [righttrigs(t)-win(1), righttrigs(t)+win(2)];
                    winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                    bins = (righttrigs(t)-win(1):binsize:righttrigs(t)+win(2));
                    binspikes = histcounts(winspikes',bins);
                    rightspikes = [rightspikes; binspikes];
                    
                end
                
                alltrigspikes{size(alltrigspikes,1)+1, 1} = [animaldayind, epoch, cell(2), cell(3)];
                alltrigspikes{size(alltrigspikes,1),2} = leftspikes; %save for later shuffling
                alltrigspikes{size(alltrigspikes,1),3} = rightspikes;
                
                %                 epochleft = [epochleft; leftspikes];
                %                 epochright = [epochright; rightspikes];
                
            end
            
        end
        
        for i = iterations:-1:1
            shuffTrialInds{animaldayind(1)}{animaldayind(2)}(:,i) = randperm(totalleft + totalright);
        end
        
        
        %         leftbinfr = (mean(epochleft,1))./binsize;
        %         rightbinfr = (mean(epochright,1))./binsize;
        %
        %         leftbinvectors = [leftbinvectors; leftbinfr];
        %         rightbinvectors = [rightbinvectors; rightbinfr];
        
        
        
    end
    
    
    % PV.(region).left = leftbinvectors;
    % PV.(region).right = rightbinvectors;
    
    
    
    %% ----- Calculate PDI -----%%
    bins =(-win(1):binsize:win(2));
    
    cells = vertcat(alltrigspikes{:,1});
    uniquecells = unique(cells(:,[1 2 4 5]),'rows');
    for c = 1:size(uniquecells,1)
        cell = uniquecells(c,:);
        
        leftspikes = vertcat(alltrigspikes{ismember(cells(:,[1 2 4 5]), cell, 'rows'),2});
        rightspikes = vertcat(alltrigspikes{ismember(cells(:,[1 2 4 5]), cell, 'rows'),3});
        
        FRleft = (mean(leftspikes,1))./binsize;
        FRleft = smoothdata(FRleft,'gaussian',5);
        FRright = (mean(rightspikes,1))./binsize;
        FRright = smoothdata(FRright,'gaussian',5);
        
        %calculate Euclidian distance. Cannot do PDI since no population
        %vectors, single cells.
        distance = (FRleft - FRright)./(FRleft + FRright);
        distance(isnan(distance))=0;
        distance = abs(distance);
        
        allspikes = vertcat(leftspikes, rightspikes);
        
        shuffDistance = nan(iterations, length(bins)-1);
        for i = 1:iterations
            if size(allspikes,1)== size(shuffTrialInds{cell(1)}{cell(2)},1)
                spikesleft_shuff = allspikes(shuffTrialInds{cell(1)}{cell(2)}(1:size(leftspikes,1),i),:);
                spikesright_shuff = allspikes(shuffTrialInds{cell(1)}{cell(2)}(size(leftspikes,1)+1:end,i),:);
            else
                shuffinds = shuffTrialInds{cell(1)}{cell(2)}(shuffTrialInds{cell(1)}{cell(2)} <= size(allspikes,1));
                shuffinds = reshape(shuffinds,size(allspikes,1),iterations);
                spikesleft_shuff = allspikes(shuffinds(1:size(leftspikes,1),i),:);
                spikesright_shuff = allspikes(shuffinds(size(leftspikes,1)+1:end,i),:);
            end
            FRleft_shuff = (mean(spikesleft_shuff,1))./binsize;
            FRleft_shuff = smoothdata(FRleft_shuff,'gaussian',5);
            FRright_shuff = (mean(spikesright_shuff,1))./binsize;
            FRright_shuff = smoothdata(FRright_shuff,'gaussian',5);
            
            shuffDistance(i,:) = (FRleft_shuff - FRright_shuff)./(FRleft_shuff + FRright_shuff);
            shuffDistance(isnan(shuffDistance))=0;
            shuffDistance = abs(shuffDistance);
        end
        
        
        for b = 1:length(bins)-1
            %                  figure
            %              histogram(shuffDistance(:,b),50);
            prc95(b) = prctile(shuffDistance(:,b),95);
            prc05(b) = prctile(shuffDistance(:,b),5);
            shuffmean(b) = mean(shuffDistance(:,b));
        end
        
        binsaftertrig = find(bins(1:end-1)>=0);
        sigbin = find((shuffmean(binsaftertrig)+prc95(binsaftertrig)) < distance(binsaftertrig), 1, 'first') + (binsaftertrig(1)-1);
        sigDistance = distance(sigbin);
        prevbin = sigbin-1;
        prevDist = distance(prevbin);
        
        [timepoint, intcpt] = polyxpoly([bins(prevbin), bins(sigbin)], [distance(prevbin), distance(sigbin)], [bins(prevbin), bins(sigbin)], [shuffmean(prevbin)+prc95(prevbin), shuffmean(prevbin)+prc95(sigbin)]);
        singleCellSelectivity.(region){c,1} = cell;
        singleCellSelectivity.(region){c,2} = timepoint;
        singleCellSelectivity.(region){c,2} = intcpt;
        
        %% Plot
        if plot == 1
            close all
            figure(r), hold on
            set(gcf, 'Position', [50,100,800,600]);
            
            plot(bins(1:end-1), distance, 'r-','LineWidth',3);
            
            %shadedErrorBar(bins(1:end-1),shuffmean,[prc95; prc05],[],1)
            patch([bins(1:end-1), fliplr(bins(1:end-1))], [shuffmean-prc05,fliplr(shuffmean)+fliplr(prc95)], 'k','FaceAlpha',0.2,'EdgeColor','none')
            p1 = plot(bins(1:end-1), shuffmean+prc95,'k-');
            p1.Color(4) = 0.25;
            p2 = plot(bins(1:end-1), shuffmean-prc05,'k-');
            p2.Color(4) = 0.25;
            plot(bins(1:end-1), shuffmean,'k-');
            
            
            ylabel('Absolute Selectivity Index')
            xlabel('Time after odor onset (seconds)');
            pause;
        end
        
        
    end
end
 save([topDir,'AnalysesAcrossAnimals\singleCellSelectivity.mat'],'singleCellSelectivity');
