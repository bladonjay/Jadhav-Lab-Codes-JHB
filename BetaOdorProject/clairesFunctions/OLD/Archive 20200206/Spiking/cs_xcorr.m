%cs_crosscorr

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir, figDir] = cs_setPaths();

bin = 0.01; % 5 ms
tmax = 0.2; % +/- 500ms for corrln
sw1 = bin*2; % for smoothing corrln.

load([topDir,'AnalysesAcrossAnimals\npCells_CA1.mat'])
allCells_CA1 = npCells;
load([topDir,'AnalysesAcrossAnimals\npCells_PFC.mat'])
allCells_PFC = npCells;


% load([topDir,'AnalysesAcrossAnimals\selectiveCells_CA1.mat'])
% allCells_CA1 = selectivecells;
% load([topDir,'AnalysesAcrossAnimals\selectiveCells_PFC.mat'])
% allCells_PFC = selectivecells;

% load([topDir,'AnalysesAcrossAnimals\plCells_CA1.mat'])
% allCells_CA1 = plcells;
% load([topDir,'AnalysesAcrossAnimals\plCells_PFC.mat'])
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
        epochs = find(~cellfun(@(x) isempty(x),nosepokeWindow{day}));
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
            
            spikesCA1 = []; spikesPFC = [];
            for ep = epochs
                npwins = nosepokeWindow{day}{ep};
                if ~isempty(spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data) && ~isempty(spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data)
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
                figure
                plot(timebase,corr_norm)
                
                %test for significance by doing a card shuffle 
                dif = (spikesCA1(end)-spikesCA1(1)).*rand(1);
                spot = spikesCA1(1)+dif;
                
                first = spikesCA1(spikesCA1 >= spot)-dif;
                last = spikesCA1(spikesCA1 < spot) - spikesCA1(1) +first(end);
                %figure, hold on
                %plot(spikesCA1,ones(1,length(spikesCA1)),'k.');

                spikesCA1 = [first;last];
                %plot(spikesCA1,ones(1,length(spikesCA1)),'r.');
                [timebase, corr_raw, corr_sm] = spiketrainxcorr(spikesCA1,spikesPFC,bin,tmax,sw1);
                corr_norm = zscore(corr_sm);
                hold on
                
                plot(timebase,corr_norm)
                close all
            end
            if any(isnan(corr_norm))
%                 disp('nan')
            else 
                allpairs = [allpairs; corr_norm];
            end
            
            
        end
    end
end

mn = mean(allpairs, 1);
[peaks,loc] = max(allpairs, [], 2);
peaktimes = timebase(loc);

meanpeak = mean(peaktimes);

[~, p] = ttest(peaktimes);
figure,
plot(timebase,mn);
hold on
plot([meanpeak, meanpeak], [-1.5 1], 'k--');
figtitle = 'CA1-PFC cross-correlation';
figfile = [figDir,'Spiking\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);
