% Replications PV differential index calculation from Igarashi et al 2012.
clear
win = [.2 1];
binsize = 0.1;


close all
%% Params
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

[topDir, figDir] = cs_setPaths();
regions = {'CA1','PFC'};

winstring = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
binstr = [num2str(binsize),'msBins'];
celltype = 'np';
switch celltype
    case 'selective'
        load([topDir, 'AnalysesAcrossAnimals\selectiveCells_CA1.mat'])
        CA1cells = selectivecells;
        load([topDir, 'AnalysesAcrossAnimals\selectiveCells_PFC.mat'])
        PFCcells = selectivecells;
    case 'np'
        load([topDir, 'AnalysesAcrossAnimals\npCells_CA1.mat'])
        CA1cells = npCells;
        load([topDir, 'AnalysesAcrossAnimals\npCells_PFC.mat'])
        PFCcells = npCells;
end

numsamples = min([size(CA1cells,1), size(PFCcells,1)]); %downsample to match region with lowest number of cells
%%
for r = 1:length(regions)
    region = regions{r};
    
    leftbinvectors = [];
    rightbinvectors = [];
    
    eval(['cells = ',region,'cells;']);
    
    %% Get big matrix with all spiking data on all trials
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
        dayleft = []; dayright = [];
        
        for ep = 1:length(runeps)
            %combine spikes over epochs
            leftspikes = [];
            rightspikes = [];
            
            epoch = runeps(ep);
            
            [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
            
            trigs = odorTriggers{day}{epoch}.allTriggers;
            
            lefttrigs = trigs(correct_left);
            righttrigs = trigs(correct_right);
            
            if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                
                for t = 1:length(lefttrigs)
                    trigwin = [lefttrigs(t)-win(1), lefttrigs(t)+win(2)+binsize];
                    winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                    bins = (lefttrigs(t)-win(1):binsize:lefttrigs(t)+win(2)+binsize);
                    binspikes = histcounts(winspikes',bins);
                    leftspikes = [leftspikes; binspikes];
                    
                end
                
                for t = 1:length(righttrigs)
                    trigwin = [righttrigs(t)-win(1), righttrigs(t)+win(2)+binsize];
                    winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
                    bins = (righttrigs(t)-win(1):binsize:righttrigs(t)+win(2)+binsize);
                    binspikes = histcounts(winspikes',bins);
                    rightspikes = [rightspikes; binspikes];
                    
                end
                
                %save cell index (animal, day, epoch, tet, cell) and spikes
                %on each trial
                
                dayleft = [dayleft; leftspikes];
                dayright = [dayright; rightspikes];
                
            else
                bins =(-win(1):binsize:win(2));
                leftspikes = nan(length(lefttrigs),length(bins));
                rightspikes = nan(length(righttrigs),length(bins));
                dayleft = [dayleft; leftspikes];
                dayright = [dayright; rightspikes];
            end
            
        end
        
        %don't save epochs
        alltrigspikes{size(alltrigspikes,1)+1, 1} = [animaldayind, cell(2), cell(3)];
        alltrigspikes{size(alltrigspikes,1),2} = dayleft; %save for later for shuffling
        alltrigspikes{size(alltrigspikes,1),3} = dayright;
        
    end
    
    %% Now that we have big matrix of spiking data, do individual bootstrap iterations
    numsamples = min([size(CA1cells,1), size(PFCcells,1)]); %downsample to match region with lowest number of cells
    
    Iterations = 500;
    crossings = [];
    for z = 1:Iterations
        
        %get random sample with replacement of cells
        disp(['Doing iteration ',num2str(z)]);
        samp = datasample(1:size(cells,1),numsamples);
        
        sampcells = alltrigspikes(samp,:);
        allcellinds = cell2mat(sampcells(:,1));
        
        %get FR PVs
        leftbinvectors= cell2mat(cellfun(@(x) nanmean(x,1)./binsize, sampcells(:,2),'UniformOutput',0));
        rightbinvectors = cell2mat(cellfun(@(x) nanmean(x,1)./binsize, sampcells(:,3),'UniformOutput',0));
        
        
        %% ----- Calculate PDI -----%%
        bins =(-win(1):binsize:win(2));
        
        PDI = zeros(1,length(bins)-1);
        for b = 1:length(bins)
            FRbinL = leftbinvectors(:,b);
            FRbinR = rightbinvectors(:,b);
            
            
            CCbin = corrcoef(FRbinL, FRbinR);
            
            PDIbin = 1-abs(CCbin(1,2)); %PV Differential Index, from Igarashi et al 2014
            PDI(b) = PDIbin;
            
        end
        
        %% ---- Shuffle -----%%
        
        iterations = 1000;
        PDI_shuff = zeros(iterations, length(bins));
        
        for i = 1:iterations
            
            %get unique animal/day inds, get number of l/r trials for each.
            %then get shuffled trials inds.
            animDays = unique(allcellinds(:,[1 2]),'rows');
            shufftrials = [];
            for j = 1:size(animDays)
                ind = animDays(j,[1 2]);
                daycellsinds = ismember(allcellinds(:,[1 2]),ind,'rows');
                daycells = sampcells(daycellsinds,:);
                numL = size(daycells{1,2},1);
                numR = size(daycells{1,3},1);
                
                trialindicator = [ones(numL,1); zeros(numR,1)];
                shuff = randsample(trialindicator, numL+numR);
                shufftrials{j,1} = shuff;
            end
            
            %now, go through all cells and shuffle trial firing rates according
            %to predefined shuffle
            for c = 1:size(sampcells,1)
                tr = [sampcells{c,2};sampcells{c,3}];
                cell = sampcells{c,1};
                %get anim/day ind
                ind = cell([1 2]);
                dayind = find(ismember(animDays,ind,'rows'));
                newL = tr(shufftrials{dayind} == 1,:);
                newR = tr(shufftrials{dayind} == 0,:);
                trigspikes_shuff{c,1} = cell;
                trigspikes_shuff{c,2} = newL;
                trigspikes_shuff{c,3} = newR;
            end
            
            %recalculate average FR and get PVs
            leftbinvectors_shuff= cell2mat(cellfun(@(x) nanmean(x,1)./binsize, trigspikes_shuff(:,2),'UniformOutput',0));
            rightbinvectors_shuff = cell2mat(cellfun(@(x) nanmean(x,1)./binsize, trigspikes_shuff(:,3),'UniformOutput',0));
            
            
            %redo PDI
            newPDI = zeros(1, length(bins));
            for b = 1:length(bins)
                FRbinL = leftbinvectors_shuff(:,b);
                FRbinR = rightbinvectors_shuff(:,b);
                
                CCbin = corrcoef(FRbinL, FRbinR);
                
                PDIbin = 1-abs(CCbin(1,2)); %PV Differential Index, from Igarashi et al 2014
                newPDI(b) = PDIbin;
                
            end
            
            PDI_shuff(i,:) = newPDI;
        end
        
        %% calculate confidence interval
        
        shuffmean = nanmean(PDI_shuff,1);
        CI = [shuffmean - prctile(PDI_shuff,5); shuffmean + prctile(PDI_shuff,95)];
        %interpolate to smooth
        newbins = bins(1):binsize/5:bins(end)-binsize/5;
        PDI = smoothdata(interp1(bins,PDI',newbins),'gaussian',15);
        shuffmean = smoothdata(interp1(bins,shuffmean',newbins),'gaussian',15);
        CI = smoothdata(interp1(bins,CI',newbins),'gaussian',15)';
        
        %% Calculate first significant bin value and save
        
        binsaftertrig = find(newbins(1:end-1)>=0);
        sigbin = find(CI(2,binsaftertrig) < PDI(binsaftertrig), 1, 'first') + (binsaftertrig(1)-1);
        
        %check to make sure PDI stays significant after this point
        while PDI(sigbin+1) < CI(2,sigbin+1)
            binstouse = binsaftertrig(find(ismember(binsaftertrig,sigbin))+1:end);
            sigbin = find((CI(2,sigbin+1)) < PDI(binstouse), 1, 'first') + (binstouse(1)-1);
        end
        sigPDI = PDI(sigbin);
        
        prevbin = sigbin-1;
        prevPDI = PDI(prevbin);
        
        %% Get time where real PDI crosses shuffle
        [timepoint, intcpt] = polyxpoly([newbins(prevbin), newbins(sigbin)], [PDI(prevbin), PDI(sigbin)], [newbins(prevbin), newbins(sigbin)], [CI(2,prevbin), CI(2,sigbin)]);
        %save distribution of timepoints
        crossings = [crossings;timepoint];
        
        
    end
    
    timepoints.(region) = crossings;
   
end
p1 = sum(timepoints.CA1 >= mean(timepoints.PFC))/length(timepoints.CA1);
    p2 = sum(timepoints.PFC >= mean(timepoints.CA1))/length(timepoints.PFC);
    p = min([p1 p2]);
    
    data = [timepoints.CA1;timepoints.PFC];
    labels = [zeros(500, 1); ones(500,1);];
    figure
    cs_boxplot(data,labels);
    
    xticklabels(regions);
    text(1, 0.6, ['p = ',num2str(p)])
    ylim([0.3 0.7])
    ylabel('PV Divergence Time (seconds)');
    figfile = [figDir,'PopSelectivity\PDI_downsample',celltype];
        %saveas(gcf,newfigfile,'fig');
        print('-dpdf', figfile);
        print('-djpeg', figfile);
    %p = sum(timepoints.CA1 >= mean(timepoints.PFC)) /length(timepoints.CA1)
    save([topDir,'AnalysesAcrossAnimals\PDI_',celltype],'timepoints')
    
