% selectivity sequences
% method 1: plot absolute value of SI, sort by peak
% method 2: plot raw FR (z score) for preferred odor, sort by peak

clear
close all

[topDir,figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
win = [0 0.8];
binsize = 0.050;
oldbins = win(1)+binsize:binsize:win(2);
newbins = win(1)+binsize:binsize/5:win(2);

for r = 1:length(regions)
    region = regions{r};
    
    %% method 1:
    %load selectivity data. already exists from cs_odorSelectivity_v2
    load([topDir,'AnalysesAcrossAnimals\selectivityData_',region])
    data = selectivityData.psth_correct; %not actually PSTH, it is selectivity index across time bins for all selective cells
    %take absolute value, find time where difference between
    %preferred/non-preferred is greatest
    data = abs(data);
    
    timebins = selectivityData.win;
    trigbins = (timebins >= 0.2 & timebins < 0.8);
    
    datatosort = data(:,trigbins);
    [~,maxinds] = max(datatosort,[],2);
    [~,cellinds] = sort(maxinds);
    data = data(cellinds,:);
    
    figure
    imagesc(timebins(trigbins),1:size(data,1),data(:,trigbins))
    figure
    imagesc(timebins,1:size(data,1),data);
    colorbar
    
    %% method 2
    
    allcells = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
        odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
        spikes = loaddatastruct(animDir, animal,'spikes');
        cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
        
        dayepochs = cs_getRunEpochs(animDir,animal,'odorplace');
        days = unique(dayepochs(:,1));
        
        for day = days'
            %% left-preferring cells
            epochs = dayepochs((dayepochs(:,1) == day),2);
            filt = ['strcmp($area,''',region,''') & strcmp($selectivity,''leftSelective'')'];
            leftcells = evaluatefilter(cellinfo{day},filt);
            leftcells = unique(leftcells(:,[2,3]),'rows');
            for c = 1:size(leftcells,1)
                cell = leftcells(c,:);
                alltrialspikes = [];
                for epoch = epochs'
                    [correct_left,~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                    trigs = odorTriggers{day}{epoch}.allTriggers(correct_left);
                    if ~isempty(spikes{day}{epoch}{cell(1)}{cell(2)}.data)
                        cellspikes = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                    else
                        continue
                    end
                    for t = 1:length(trigs)
                        trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
                        winspikes = cellspikes(isExcluded(cellspikes,trigwin));
                        bins = (trigs(t)-win(1):binsize:trigs(t)+win(2));
                        binspikecount = histcounts(winspikes,bins);
                        alltrialspikes = [alltrialspikes;binspikecount];
                    end
                end
                binfr = mean(alltrialspikes,1)./binsize;
                zscorebinfr = (binfr-mean(binfr))/std(binfr);
                interpbinfr = interp1(oldbins,zscorebinfr,newbins);
                allcells = [allcells;interpbinfr];
            end
            
            %% right-preferring cells
            filt = ['strcmp($area,''',region,''') & strcmp($selectivity,''rightSelective'')'];
            rightcells = evaluatefilter(cellinfo{day},filt);
            rightcells = unique(rightcells(:,[2,3]),'rows');
            for c = 1:size(rightcells,1)
                cell = rightcells(c,:);
                alltrialspikes = [];
                for epoch = epochs'
                    [~,correct_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                    trigs = odorTriggers{day}{epoch}.allTriggers(correct_right);
                    if ~isempty(spikes{day}{epoch}{cell(1)}{cell(2)}.data)
                        cellspikes = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                    else
                        continue
                    end
                    for t = 1:length(trigs)
                        trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
                        winspikes = cellspikes(isExcluded(cellspikes,trigwin));
                        bins = (trigs(t)-win(1):binsize:trigs(t)+win(2));
                        binspikecount = histcounts(winspikes,bins);
                        alltrialspikes = [alltrialspikes;binspikecount];
                    end
                end
                binfr = mean(alltrialspikes,1)./binsize;
                zscorebinfr = (binfr-mean(binfr))/std(binfr);
                interpbinfr = interp1(oldbins,zscorebinfr,newbins);
                allcells = [allcells;interpbinfr];
            end
        end
    end
    
    trigbins = (newbins >= 0 & newbins < 0.8);
    
    datatosort = allcells(:,trigbins);
    [~,maxinds] = max(datatosort,[],2);
    [s,cellinds] = sort(maxinds);
    data = allcells(cellinds,:);
    
    imagesc(newbins,1:size(data,1),data)
    colorbar
    
    xlabel('Time from odor onset (s)');
    ylabel('Cell number');
    title([region, ' selective cells firing rate'])
    
    
    figfile = [figDir,'OdorSelectivity\',region,' selective cell firing rate'];
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
    
end