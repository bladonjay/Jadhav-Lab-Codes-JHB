%cs_xcov
%use DFAsj_getthetacrosscov


clear
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

load([topDir,'AnalysesAcrossAnimals\npCells_CA1.mat'])
allCellsCA1 = npCells;
load([topDir,'AnalysesAcrossAnimals\npCells_PFC.mat'])
allCellsPFC = npCells;
filetag = 'np';

% 
% load([topDir,'AnalysesAcrossAnimals\selectiveCells_CA1.mat'])
% allCellsCA1 = selectivecells;
% load([topDir,'AnalysesAcrossAnimals\selectiveCells_PFC.mat'])
% allCellsPFC = selectivecells;
% filetag = 'selective';

%% -- All Cells
    cov_allpairs = [];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal,'Expt\',animal,'_direct\'];
        
        %load data
        npWins = loaddatastruct(animDir,animal,'nosepokeWindow');
        cellinfo = loaddatastruct(animDir,animal,'cellinfo');
        spikes = loaddatastruct(animDir, animal, 'spikes');

        %get days and epochs
        runeps = cs_getRunEpochs(animDir, animal,'odorplace');
        days = unique(runeps(:,1));
        
        for day = days'
            
            epochs = cs_getRunEpochs(animDir, animal,'odorplace',day);
            epochs = epochs(:,2);
            %eps = runeps(runeps(:,1) == day,2);
            
            %get cell pairs
            mat = [a,day];
            cellsCA1 = allCellsCA1(ismember(allCellsCA1(:,[1,2]),mat,'rows'),[3,4]);
            
            cellsPFC = allCellsPFC(ismember(allCellsPFC(:,[1,2]),mat,'rows'),[3,4]);
            
            if isempty(cellsCA1) || isempty(cellsPFC)
                continue
            end
            disp(['Doing ',animal,' day ',num2str(day)]);
            [A,B] = meshgrid(1:size(cellsCA1,1),1:size(cellsPFC,1));
            pairInds=reshape(cat(2,A',B'),[],2);
            
            for p = 1:size(pairInds,1)
                c1 = cellsCA1(pairInds(p,1),:);
                c2 = cellsPFC(pairInds(p,2),:);
                
                eps1 = cs_findGoodEpochs(spikes{day},{'data'},c1);
                eps2 = cs_findGoodEpochs(spikes{day},{'data'},c2);
                eps = intersect(epochs,intersect(eps1,eps2));
                
                cov = [];
                for ep = eps'
                    ind = [day, ep, cellsCA1(pairInds(p,1),:), cellsPFC(pairInds(p,2),:)];
                    
                    %find full time range
                    timerange = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange;
                    samprate = 1500;
                    times = timerange(1):1/samprate:timerange(2);
                    totaltime = diff(timerange);
                    
                    
                    %get NP windows
                    wins = npWins{day}{ep};
                    
                    %must get 'exclude periods' for cross cov function - i.e. any times NOT in NP, function will then use NP periods only
                    [wintime,winvec] = wb_list2vec(wins,times);
                    nonnpvec = ~winvec;
                    nonnplist = vec2list(nonnpvec,times);
                    

                    %get cross covariance for each cell pair
                    out = DFAsj_getthetacrosscov(ind, nonnplist, spikes);
                    
                    if any(~isnan(out.Zcrosscov)) && out.Neventscorr > 50
                        timebins = out.corr.time(out.bins);
                        cov = [cov; out.Zcrosscov];
                    end
                end
                 mncov = mean(cov,1);
                cov_allpairs = [cov_allpairs;mncov];
                
                
            end
        end
    end
    
   mn = mean(cov_allpairs, 1);
mn_sm = smoothdata(mn,'gaussian',8);
sem = stderr(cov_allpairs);
sem_sm = smoothdata(sem,'gaussian',8);
figure,

timebase = out.corr.time;
patch([timebase, fliplr(timebase)],[mn_sm+sem_sm, fliplr(mn_sm-sem_sm)],'b') 
hold on
plot(timebase, mn_sm,'k')

[peaks,loc] = max(cov_allpairs, [], 2);
peaktimes = timebase(loc);

meanpeak = mean(peaktimes);

[~, p] = ttest(peaktimes);
% figure,
% plot(timebase,mn);
% hold on
plot([meanpeak, meanpeak], [min(mn) max(mn)], 'k--');
figtitle = ['CA1-PFC xcov ',filetag];
figfile = [figDir,'Spiking\',figtitle];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);

save([topDir,'AnalysesAcrossAnimals\xcov_',filetag],'cov_allpairs')