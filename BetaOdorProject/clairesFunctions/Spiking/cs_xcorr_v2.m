%cs_crosscorr
%makes histogram of peak times
close all
clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir, figDir] = cs_setPaths();

bin = 0.01; % 10 ms
tmax = 0.2; % +/- 500ms for corrln
sw1 = bin*2; % for smoothing corrln.

load([topDir,'AnalysesAcrossAnimals\npCells_CA1.mat'])
allCells_CA1 = npCells;
load([topDir,'AnalysesAcrossAnimals\npCells_PFC.mat'])
allCells_PFC = npCells;
% 

% load([topDir,'AnalysesAcrossAnimals\selectiveCells_CA1.mat'])
% allCells_CA1 = selectivecells;
% load([topDir,'AnalysesAcrossAnimals\selectiveCells_PFC.mat'])
% allCells_PFC = selectivecells;

% load([topDir,'AnalysesAcrossAnimals\plCells_CA1-OB.mat'])
% allCells_CA1 = plcells;
% load([topDir,'AnalysesAcrossAnimals\plCells_PFC-OB.mat'])
% allCells_PFC = plcells;

animdaysCA1 = unique(allCells_CA1(:,[1 2]),'rows');
animdaysPFC = unique(allCells_PFC(:,[1 2]),'rows');

animdays = intersect(animdaysCA1, animdaysPFC,'rows');
animinds = unique(animdays(:,1));

allpairs = [];
for a = 1:length(animinds)
    animind = animinds(a);
    animal = animals{animind};
    animDir = [topDir,animal,'Expt\', animal,'_direct\'];
    days = animdays(animdays(:,1) == animind, 2);
    
    for d = 1:length(days)
        day = days(d);
        daystring = getTwoDigitNumber(day);
        load([animDir, animal,'spikes',daystring,'.mat']);
        load([animDir, animal,'nosepokeWindow',daystring,'.mat']);
        
        epochs = cs_getRunEpochs(animDir, animal,'odorplace',day);
        epochs = epochs(:,2);
        
        disp(['Doing animal ', animal, ' day ', daystring]);
        
        %cells from the day
        cellsCA1 = allCells_CA1(ismember(allCells_CA1(:,[1,2]),[animind, day],'rows'),[3,4]);
        cellsPFC = allCells_PFC(ismember(allCells_PFC(:,[1,2]),[animind, day],'rows'),[3,4]);
        
        %create pairs of cells
        [A,B] = meshgrid(1:size(cellsCA1,1),1:size(cellsPFC,1));
        pairInds=reshape(cat(2,A',B'),[],2);
        
        for p = 1:size(pairInds,1)
            cCA1 = cellsCA1(pairInds(p,1),:);
            cPFC = cellsPFC(pairInds(p,2),:);
            
            eps1 = cs_findGoodEpochs(spikes{day},{'data'},cCA1);
            eps2 = cs_findGoodEpochs(spikes{day},{'data'},cPFC);
            eps = intersect(epochs,intersect(eps1,eps2));
            
            spikesCA1 = []; spikesPFC = [];
            for ep = eps'
                npwins = nosepokeWindow{day}{ep};
                if ~isempty(spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data )&& ~isempty(spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data)
                    rawspikesCA1 = spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data(:,1);
                    rawspikesPFC = spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data(:,1);
                    
                    
                    %find spikes that occured during nosepoke
                    spikesCA1 = [spikesCA1; rawspikesCA1(any(rawspikesCA1' >= npwins(:,1) & rawspikesCA1' < npwins(:,2),1))];
                    spikesPFC = [spikesPFC; rawspikesPFC(any(rawspikesPFC' >= npwins(:,1) & rawspikesPFC' < npwins(:,2),1))];
                end
            end
            
            if ~isempty(spikesCA1) && ~isempty(spikesPFC)
                [timebase, corr_raw, corr_sm] = spiketrainxcorr(spikesCA1,spikesPFC,bin,tmax,sw1);
                corr_norm = zscore(corr_sm);
            end
            if any(isnan(corr_norm))
%                 disp('nan')
            else 
                
                
                binsize = timebase(2)-timebase(1);
                newbins = timebase(1):binsize/5:timebase(end)-binsize/5;
                newcorr = smoothdata(interp1(timebase,corr_norm',newbins),'gaussian',10)';
                figure,plot(newbins,newcorr);
                hold on
                plot([0 0],[0 max(newcorr)],'k--');
                close
                
                [~,peak] = max(corr_norm);
                peaktime = timebase(peak);
                
                allpairs = [allpairs; peaktime];
            end
            
            
        end
    end
end
mn = mean(allpairs);

[counts,edges] = histcounts(allpairs,20,'Normalization','probability');
binsize = edges(2)-edges(1);
bins =edges(2:end)-(binsize/2);

newbins = bins(1):binsize/5:bins(end)-binsize/5;

newcounts = smoothdata(interp1(bins,counts',newbins),'gaussian',10)';

figure,
plot(newbins,newcounts);
hold on
plot([0 0],[0 0.14],'k--');
% histogram(allpairs, 40);
% hold on
% plot(mn, 0, 'r*', 'MarkerSize',10)
figtitle = 'CA1-PFC cross-correlation all';
figfile = [figDir,'Spiking\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);


%% ------ Selective Cells

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir, figDir] = cs_setPaths();

bin = 0.01; % 10 ms
tmax = 0.2; % +/- 500ms for corrln
sw1 = bin*2; % for smoothing corrln.


load([topDir,'AnalysesAcrossAnimals\selectiveCells_CA1.mat'])
allCells_CA1 = selectivecells;
load([topDir,'AnalysesAcrossAnimals\selectiveCells_PFC.mat'])
allCells_PFC = selectivecells;

animdaysCA1 = unique(allCells_CA1(:,[1 2]),'rows');
animdaysPFC = unique(allCells_PFC(:,[1 2]),'rows');

animdays = intersect(animdaysCA1, animdaysPFC,'rows');
animinds = unique(animdays(:,1));

allpairs = [];
for a = 1:length(animinds)
    animind = animinds(a);
    animal = animals{animind};
    animDir = [topDir,animal,'Expt\', animal,'_direct\'];
    days = animdays(animdays(:,1) == animind, 2);
    
    for d = 1:length(days)
        day = days(d);
        daystring = getTwoDigitNumber(day);
        load([animDir, animal,'spikes',daystring,'.mat']);
        load([animDir, animal,'nosepokeWindow',daystring,'.mat']);
         epochs = cs_getRunEpochs(animDir, animal,'odorplace',day);
        epochs = epochs(:,2);
        
        disp(['Doing animal ', animal, ' day ', daystring]);
        
        %cells from the day
        cellsCA1 = allCells_CA1(ismember(allCells_CA1(:,[1,2]),[animind, day],'rows'),[3,4]);
        cellsPFC = allCells_PFC(ismember(allCells_PFC(:,[1,2]),[animind, day],'rows'),[3,4]);
        
        %create pairs of cells
        [A,B] = meshgrid(1:size(cellsCA1,1),1:size(cellsPFC,1));
        pairInds=reshape(cat(2,A',B'),[],2);
        
        for p = 1:size(pairInds,1)
            cCA1 = cellsCA1(pairInds(p,1),:);
            cPFC = cellsPFC(pairInds(p,2),:);
            
             eps1 = cs_findGoodEpochs(spikes{day},{'data'},cCA1);
            eps2 = cs_findGoodEpochs(spikes{day},{'data'},cPFC);
            eps = intersect(epochs,intersect(eps1,eps2));
            
            spikesCA1 = []; spikesPFC = [];
            for ep = eps'
                npwins = nosepokeWindow{day}{ep};
                if ~isempty(spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data )&& ~isempty(spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data)
                    rawspikesCA1 = spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data(:,1);
                    rawspikesPFC = spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data(:,1);
                    
                    
                    %find spikes that occured during nosepoke
                    spikesCA1 = [spikesCA1; rawspikesCA1(any(rawspikesCA1' >= npwins(:,1) & rawspikesCA1' < npwins(:,2),1))];
                    spikesPFC = [spikesPFC; rawspikesPFC(any(rawspikesPFC' >= npwins(:,1) & rawspikesPFC' < npwins(:,2),1))];
                end
            end
            
            if ~isempty(spikesCA1) && ~isempty(spikesPFC)
                [timebase, corr_raw, corr_sm] = spiketrainxcorr(spikesCA1,spikesPFC,bin,tmax,sw1);
                corr_norm = zscore(corr_sm);
            end
            if any(isnan(corr_norm))
%                 disp('nan')
            else 
                [~,peak] = max(corr_norm);
                peaktime = timebase(peak);
                
                allpairs = [allpairs; peaktime];
            end
            
            
        end
    end
end
mn = mean(allpairs);

[counts,edges] = histcounts(allpairs,20,'Normalization','probability');
binsize = edges(2)-edges(1);
bins =edges(2:end)-(binsize/2);

newbins = bins(1):binsize/5:bins(end)-binsize/5;

newcounts = smoothdata(interp1(bins,counts',newbins),'gaussian',10)';

figure,
plot(newbins,newcounts);
hold on
plot([0 0],[0 0.14],'k--');
% histogram(allpairs, 40);
% hold on
% plot(mn, 0, 'r*', 'MarkerSize',10)
figtitle = 'CA1-PFC cross-correlation selective';
figfile = [figDir,'Spiking\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);

