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

%%
for r = 1:length(regions)
    region = regions{r};
    
    leftbinvectors = [];
    rightbinvectors = [];
    
%     load selective cells only
    load([topDir, 'AnalysesAcrossAnimals\selectiveINs_', region, '.mat'])
    cells = selectiveint;

%     load([topDir, 'AnalysesAcrossAnimals\npCells_', region, '.mat'])
%     cells = npCells;
%         
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
        %combine spikes over epochs
        for ep = 1:length(runeps)
            epoch = runeps(ep);
            
            leftspikes = [];
            rightspikes = [];
            if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);

                [correct_left, correct_right, ~, ~] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                
                trigs = odorTriggers{day}{epoch}.allTriggers;

                lefttrigs = trigs(correct_left);
                righttrigs = trigs(correct_right);

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
                
                alltrigspikes{size(alltrigspikes,1)+1, 1} = [animaldayind, epoch, cell(2), cell(3)];
                alltrigspikes{size(alltrigspikes,1),2} = leftspikes; %save for later for shuffling
                alltrigspikes{size(alltrigspikes,1),3} = rightspikes;
             
                epochleft = [epochleft; leftspikes];
                epochright = [epochright; rightspikes];
                
            end
            


        end
        
        
        leftbinfr = (mean(epochleft,1))./binsize;
        %leftbinfr = smoothdata(leftbinfr);
        
        rightbinfr = (mean(epochright,1))./binsize;
        %rightbinfr = smoothdata(rightbinfr);

        leftbinvectors = [leftbinvectors; leftbinfr];
        rightbinvectors = [rightbinvectors; rightbinfr];
       
             
        
    end


PV.(region).left = leftbinvectors;
PV.(region).right = rightbinvectors;


        
%% ----- Calculate PDI -----%%
bins =(-win(1):binsize:win(2));
                       
PDI = zeros(1,length(bins));
for b = 1:length(bins)
    FRbinL = leftbinvectors(:,b);
    FRbinR = rightbinvectors(:,b);
    

    CCbin = corrcoef(FRbinL, FRbinR);

    PDIbin = 1-abs(CCbin(1,2)); %PV Differential Index, from Igarashi et al 2014
    PDI(b) = PDIbin;
    
end



%% ---- Shuffle -----%%

iterations = 1000;
allcellinds = {alltrigspikes{:,1}}; allcellinds = allcellinds';
allcellinds = cell2mat(allcellinds);

sizeL = cell2mat(cellfun(@(x) size(x,1),alltrigspikes(:,2),'UniformOutput',false));
sizeR = cell2mat(cellfun(@(x) size(x,1),alltrigspikes(:,3),'UniformOutput',false));

metainds = [allcellinds(:,[1 2 3]), sizeL, sizeR];

metainds = unique(metainds,'rows');
% cols 1-3 are animal,day,epoch index; col 4 is number of left trials, col
% 5 is number of right trials. 
alltrigspikesind = cell2mat({alltrigspikes{:,1}}');
alltrigspikesind = alltrigspikesind(:,[1 2 4 5]);


shuffPDI = zeros(iterations, length(bins));

for i = 1:iterations
    if rem(i,100) == 0
    disp (['Doing shuffle iteration ',num2str(i)]);
    end
    for t = 1:size(metainds,1)
        
        lv = ones(metainds(t,4),1);
        rv = zeros(metainds(t,5),1);
        
        %use randsample since we are just shuffling the order
        all = [lv; rv];
        shufftrigs = randsample(all, size(all,1));
        allshufftrigs{t,1} = shufftrigs;  
        
    end
    
    newleftbinvectors = [];
    newrightbinvectors = [];
    
    for c = 1:size(cells,1)
        cell = cells(c,:);
        animaldayind = cells(c,[1 2]);
        epochs = allcellinds(ismember(allcellinds(:,[1 2 4 5]), cell, 'rows'),3);
        
        epochinds = find(ismember(metainds(:,[1 2 3]),[repmat(animaldayind, length(epochs),1), epochs],'rows'));
        
        
        cellinds = find(ismember(alltrigspikesind, cell, 'rows')); %number of elements should be the same as number of epochs
            if ~isequal(length(cellinds), length(epochs))
                n = cell2mat(alltrigspikes(cellinds));
                epochs = n(:,3);
            end
                
        newleft = [];
        newright = [];
        for e = 1:length(epochs)
            epochtrigs = allshufftrigs{epochinds(e)};
            newLinds = find(epochtrigs == 1);
            newRinds = find(epochtrigs == 0);
            
            leftspikes = alltrigspikes{cellinds(e), 2};
            rightspikes = alltrigspikes{cellinds(e), 3};
            allspikes = [leftspikes; rightspikes];
            
            newleft = [newleft; allspikes(newLinds, :)];
            newright = [newright; allspikes(newRinds, :)];
            
        end
        
        newleftbinfr = (mean(newleft,1))./binsize;
        newrightbinfr = (mean(newright,1))./binsize;

        newleftbinvectors = [newleftbinvectors; newleftbinfr];
        newrightbinvectors = [newrightbinvectors; newrightbinfr];
    end
    
    
 %% --- calculate new PDI
                     
    newPDI = zeros(1,length(bins));
    for b = 1:length(bins)
        FRbinL = newleftbinvectors(:,b);
        FRbinR = newrightbinvectors(:,b);

%         magL = sqrt(sum(FRbinL.^2));
%         magR = sqrt(sum(FRbinR.^2));
% 
%         CCbin = (dot(FRbinL,FRbinR))/(magL * magR);

        CCbin = corrcoef(FRbinL, FRbinR);

        PDIbin = 1-CCbin(1,2); %PV Differential Index, from Igarashi et al 2014
        newPDI(b) = PDIbin;

    end
    
    shuffPDI(i,:) = newPDI;

end

%% calculate confidence interval
    shuffmean = mean(shuffPDI,1);
    CI = [shuffmean - prctile(shuffPDI,5); shuffmean + prctile(shuffPDI,95)];
    %interpolate to smooth
    newbins = bins(1):binsize/5:bins(end)-binsize/5;
    PDI = smoothdata(interp1(bins,PDI',newbins),'gaussian',15);
    shuffmean = smoothdata(interp1(bins,shuffmean',newbins),'gaussian',15);
    CI = smoothdata(interp1(bins,CI',newbins),'gaussian',15)';
    %CI = getCI(shuffPDI,95);

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
    
    [timepoint, intcpt] = polyxpoly([newbins(prevbin), newbins(sigbin)], [PDI(prevbin), PDI(sigbin)], [newbins(prevbin), newbins(sigbin)], [CI(2,prevbin), CI(2,sigbin)]); 
    significantPDI.(region) = [timepoint, intcpt];

%% Plot
    figure(r), hold on
    set(gcf, 'Position', [50,100,800,600]);
    
    %PDI = smoothdata(PDI);
    plot(newbins, PDI, 'k-','LineWidth',3);
    
    maximums(r) = max(PDI);

    %shadedErrorBar(bins(1:end-1),shuffmean,[prc95; prc05],[],1)
    patch([newbins, fliplr(newbins)], [CI(1,:),fliplr(CI(2,:))], 'k','FaceAlpha',0.2,'EdgeColor','none')
    %p1 = plot(bins, shuffmean+prc95,'k-');
    %p1.Color(4) = 0.25;
    %p2 = plot(bins, shuffmean-prc05,'k-');
    %p2.Color(4) = 0.25;
    %plot(bins, shuffmean,'k-');
    
    text(0.2, 0.7, ['time = ', num2str(round(timepoint,2)),' s']);
    
    
    
    ylabel('PDI')
    xlabel('Time after odor onset (seconds)');
    title(region);
    
end
save([topDir, 'AnalysesAcrossAnimals\significantPDI_int.mat'], 'significantPDI');
save([topDir, 'AnalysesAcrossAnimals\PV_int.mat'], 'PV');

for r = 1:length(regions)
    region = regions{r};
    figure(r)
    axis([-win(1) win(2) 0 max(maximums)+0.05]);
    plot([0 0], [0 max(maximums)+0.05], 'k--', 'LineWidth',3);
    
    figtitle = ['EnsemblePDI_IN',region];
    
    figfile = [figDir,'Interneurons\',figtitle];
    
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    %saveas(gcf,figfile,'fig');
    %legend({region, 'Shuffle'})
end