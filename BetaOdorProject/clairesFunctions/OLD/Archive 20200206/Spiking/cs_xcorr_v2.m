%cs_crosscorr
%Take an average of the entire trace- don't find peaks
clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir, figDir] = cs_setPaths();

bin = 0.01; % 10 ms
tmax = 0.1; % +/- 500ms for corrln
sw1 = bin*2; % for smoothing corrln.

celltype = 'choice';

switch celltype
    case 'all'
        load([topDir,'AnalysesAcrossAnimals\npCells_CA1.mat'])
        allCells_CA1 = npCells;
        load([topDir,'AnalysesAcrossAnimals\npCells_PFC.mat'])
        allCells_PFC = npCells;
        
    case 'selective'
        load([topDir,'AnalysesAcrossAnimals\selectiveCells_CA1.mat'])
        allCells_CA1 = selectivecells;
        load([topDir,'AnalysesAcrossAnimals\selectiveCells_PFC.mat'])
        allCells_PFC = selectivecells;
        
    case 'phaselocked'
%         load([topDir,'AnalysesAcrossAnimals\PhaseLocking\plCells_beta_CA1-OB.mat'])
%         allCells_CA1 = plcells;
        load([topDir,'AnalysesAcrossAnimals\PhaseLocking\plCells_beta_CA1-CA1.mat'])
        allCells_CA1 = plcells;
        load([topDir,'AnalysesAcrossAnimals\PhaseLocking\plCells_beta_CA1-CA1.mat'])
        allCells_CA1 = [allCells_CA1;plcells];
        allCells_CA1 = unique(allCells_CA1,'rows');
        
        load([topDir,'AnalysesAcrossAnimals\PhaseLocking\plCells_beta_PFC-CA1.mat'])
        allCells_PFC = plcells;
        load([topDir,'AnalysesAcrossAnimals\PhaseLocking\plCells_beta_PFC-PFC.mat'])
        allCells_PFC = [allCells_PFC; plcells];
        allCells_PFC = unique(allCells_PFC,'rows');
    case 'odor'
        load([topDir,'AnalysesAcrossAnimals\odorCells_CA1.mat'])
        allCells_CA1 = odorcells;
        load([topDir,'AnalysesAcrossAnimals\odorCells_PFC.mat'])
        allCells_PFC = odorcells;
        
    case 'choice'
        load([topDir,'AnalysesAcrossAnimals\choiceCells_CA1.mat'])
        allCells_CA1 = choicecells;
        load([topDir,'AnalysesAcrossAnimals\choiceCells_PFC.mat'])
        allCells_PFC = choicecells;
        
end

animdaysCA1 = unique(allCells_CA1(:,[1 2]),'rows');
animdaysPFC = unique(allCells_PFC(:,[1 2]),'rows');

animdays = intersect(animdaysCA1, animdaysPFC,'rows');
animinds = unique(animdays(:,1));

allpairs = [];
for a = 1:length(animinds)
    animind = animinds(a);
    animal = animals{animind};
    animDir = [topDir,animal,'Expt\', animal,'_direct\'];
    load([animDir, animal,'highBeta.mat']);
    days = animdays(animdays(:,1) == animind, 2);
    
    for d = 1:length(days)
        day = days(d);
        daystring = getTwoDigitNumber(day);
        load([animDir, animal,'spikes',daystring,'.mat']);
        %load([animDir, animal,'nosepokeWindow',daystring,'.mat']);
        epochs = find(~cellfun(@(x) isempty(x),highBeta{day}));
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
                %npwins = nosepokeWindow{day}{ep};
                betawins = highBeta{day}{ep}.OB;
                if ~isempty(spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data) && ~isempty(spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data)
                    rawspikesCA1 = spikes{day}{ep}{cCA1(1)}{cCA1(2)}.data(:,1);
                    rawspikesPFC = spikes{day}{ep}{cPFC(1)}{cPFC(2)}.data(:,1);
                    
                    for tr = 1:size(betawins,1)
                        spikesCA1 = rawspikesCA1(isExcluded(rawspikesCA1,betawins(tr,:)));
                        spikesPFC = rawspikesPFC(isExcluded(rawspikesPFC,betawins(tr,:)));
                        if ~isempty(spikesCA1) && ~isempty(spikesPFC)
                            %find spikes that occured during high beta
                            %                     spikesCA1 = [spikesCA1; rawspikesCA1(any(rawspikesCA1' >= betawins(:,1) & rawspikesCA1' < betawins(:,2),1))];
                            %                     spikesPFC = [spikesPFC; rawspikesPFC(any(rawspikesPFC' >= betawins(:,1) & rawspikesPFC' < betawins(:,2),1))];
                            [timebase, corr_raw, corr_sm] = spiketrainxcorr(spikesPFC,spikesCA1,bin,tmax,sw1);
                            corr_norm = zscore(corr_sm);
                            if ~isnan(corr_norm)
                            allpairs = [allpairs; corr_norm];
                            end
                        end
                    end
                end
            end
            
        end
        
    end
    
end

mn = mean(allpairs(:,2:end-1),1);
sem = stderr(allpairs(:,2:end-1));
figure,
%histogram(allpairs, 40);
plot(timebase(2:end-1),mn)
hold on
patch([timebase(2:end-1), fliplr(timebase(2:end-1))], [mn-sem,fliplr(mn)+fliplr(sem)], 'k','FaceAlpha',0.2,'EdgeColor','none')

plot([0 0], [-.15 0.35], 'k--')
figtitle = ['CA1-PFC cross-correlation ', celltype];
figfile = [figDir,'Spiking\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);
