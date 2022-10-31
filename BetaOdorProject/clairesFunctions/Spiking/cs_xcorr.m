%cs_crosscorr
close all

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir, figDir] = cs_setPaths();

bin = 0.01; % 10 ms
tmax = 0.25; % +/- 500ms for corrln
sw1 = bin*2; % for smoothing corrln.

% load([topDir,'AnalysesAcrossAnimals\npCells_CA1.mat'])
% allCells_CA1 = npCells;
% load([topDir,'AnalysesAcrossAnimals\npCells_PFC.mat'])
% allCells_PFC = npCells;
% filetag = 'np';


load([topDir,'AnalysesAcrossAnimals\selectiveCells_CA1.mat'])
allCells_CA1 = selectivecells;
load([topDir,'AnalysesAcrossAnimals\selectiveCells_PFC.mat'])
allCells_PFC = selectivecells;
filetag = 'selective';

% load([topDir,'AnalysesAcrossAnimals\plCells_CA1.mat'])
% allCells_CA1 = plcells;
% load([topDir,'AnalysesAcrossAnimals\plCells_PFC.mat'])
% allCells_PFC = plcells;

animdaysCA1 = unique(allCells_CA1(:,[1 2]),'rows');
animdaysPFC = unique(allCells_PFC(:,[1 2]),'rows');

animdays = intersect(animdaysCA1, animdaysPFC,'rows');
animinds = unique(animdays(:,1));

allpairs = [];
allpairs_raw = [];

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
                if ~isempty(spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data) && ~isempty(spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data)
                    rawspikesCA1 = spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data(:,1);
                    rawspikesPFC = spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data(:,1);
                    
                    
                    %find spikes that occured during nosepoke
                     spikesCA1 = [spikesCA1; rawspikesCA1(isExcluded(rawspikesCA1, npwins))];
                     spikesPFC = [spikesPFC ;rawspikesPFC(isExcluded(rawspikesPFC, npwins))];
                end
            end
            
            
            if ~isempty(spikesCA1) && ~isempty(spikesPFC)
                [timebase, corr_raw, corr_sm] = spiketrainxcorr(spikesCA1,spikesPFC,bin,tmax,sw1);
                corr_norm = corr_sm;
                corr = corr_raw;
            end
            if any(isnan(corr_norm))
%                 disp('nan')
            else 
                allpairs = [allpairs; corr_norm];
                allpairs_raw = [allpairs_raw; corr_raw];
            end
            
            
        end
    end
end

mn = mean(allpairs_raw, 1);
mn_sm = smoothdata(mn,'gaussian',8);
sem = stderr(allpairs_raw);
sem_sm = smoothdata(sem,'gaussian',8);
figure,

patch([timebase, fliplr(timebase)],[mn_sm+sem_sm, fliplr(mn_sm-sem_sm)],'b') 
hold on
plot(timebase, mn_sm,'k')

[peaks,loc] = max(allpairs, [], 2);
peaktimes = timebase(loc);

meanpeak = mean(peaktimes);

[~, p] = ttest(peaktimes);
% figure,
% plot(timebase,mn);
% hold on
plot([0, 0], [min(mn) max(mn)], 'k--');
figtitle =['CA1-PFC xcorr ',filetag];
figfile = [figDir,'Spiking\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);

save([topDir,'AnalysesAcrossAnimals\xcorr_',filetag],'allpairs_raw','timebase')
clear
close all