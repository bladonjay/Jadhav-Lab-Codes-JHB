%cs_xcov
%use DFAsj_getthetacrosscov


clear
[topDir, figDir] = cs_setPaths;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

allcells = 0;
selectivecells = 0;
sameprefcells = 1;

%% -- All Cells
if allcells == 1
cov = [];
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
        
        eps = runeps(runeps(:,1) == day,2);
        
        for ep = eps'
            
            %find first good tet
            filt = '$numspikes > 0';
            tets = evaluatefilter(cellinfo{day}{ep},filt);
            tet = tets(1,1);
            
            eeg = loadeegstruct(animDir, animal, 'eeg',day,ep,tet);
            times = geteegtimes(eeg{day}{ep}{tet});
            
            %get NP windows
            wins = npWins{day}{ep};
            
            %must get 'exclude periods' for cross cov function - i.e. any times NOT in NP, function will then use NP periods only
            [wintime,winvec] = wb_list2vec(wins,times);
            nonnpvec = ~winvec;
            nonnplist = vec2list(nonnpvec,times);
            
            %get cell pairs
            filt = 'strcmp($area,''CA1'')';
            cellsCA1 = evaluatefilter(cellinfo{day}{ep},filt);
            
            filt = 'strcmp($area,''PFC'')';
            cellsPFC = evaluatefilter(cellinfo{day}{ep},filt);
            
            [A,B] = meshgrid(1:size(cellsCA1,1),1:size(cellsPFC,1));
            pairInds=reshape(cat(2,A',B'),[],2);
            
            for p = 1:size(pairInds,1)
                
                ind = [day, ep, cellsCA1(pairInds(p,1),:), cellsPFC(pairInds(p,2),:)];
                
                %get cross covariance for each cell pair
                out = DFAsj_getthetacrosscov(ind, nonnplist, spikes);
                
                if any(~isnan(out.Zcrosscov_sm)) && out.Neventscorr > 10
                    timebins = out.corr.time(out.bins);
                    cov = [cov; out.Zcrosscov_sm];
                end
            end

            
        end
    end
end
[~,maxcovind] = max(cov');
%t = -0.2:0.01:0.2;
%timebins = out.corr.time(out.bins);
maxcov= timebins(maxcovind);



t = -0.2:0.02:0.2;
%h = histcounts(maxcov,t,'Normalization','probability');
h = histcounts(maxcov,t);
figure
bar(t(2:end),h);

hold on
plot([0 0],[0 max(h)],'k--')

xlabel('Peak cross-covariance time (seconds)')
ylabel('Fraction of cell pairs');

figFile = [figDir,'Spiking\crosscov_CA1-PFC_allcells'];
print(figFile,'-djpeg');
print(figFile,'-dpdf');
end
%% --- Selective Cells Only
if selectivecells == 1
cov = [];
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
        
        eps = runeps(runeps(:,1) == day,2);
        
        for ep = eps'
            
            %find first good tet
            filt = '$numspikes > 0';
            tets = evaluatefilter(cellinfo{day}{ep},filt);
            tet = tets(1,1);
            
            eeg = loadeegstruct(animDir, animal, 'eeg',day,ep,tet);
            times = geteegtimes(eeg{day}{ep}{tet});
            
            
            wins = npWins{day}{ep};
            
            [wintime,winvec] = wb_list2vec(wins,times);
            nonnpvec = ~winvec;
            nonnplist = vec2list(nonnpvec,times);
            
            
            filt = 'strcmp($area,''CA1'') && (strcmp($selectivity,''leftSelective'') || strcmp($selectivity,''rightSelective''))';
            cellsCA1 = evaluatefilter(cellinfo{day}{ep},filt);
            
            
            filt = 'strcmp($area,''PFC'') && (strcmp($selectivity,''leftSelective'') || strcmp($selectivity,''rightSelective''))';
            cellsPFC = evaluatefilter(cellinfo{day}{ep},filt);
            
            [A,B] = meshgrid(1:size(cellsCA1,1),1:size(cellsPFC,1));
            pairInds=reshape(cat(2,A',B'),[],2);
            
            for p = 1:size(pairInds,1)
                
                ind = [day, ep, cellsCA1(pairInds(p,1),:), cellsPFC(pairInds(p,2),:)];
                
                out = DFAsj_getthetacrosscov(ind, nonnplist, spikes);
                
                if any(~isnan(out.Zcrosscov_sm)) && out.Neventscorr > 10
                cov = [cov; out.Zcrosscov_sm];
                end
            end

            
        end
    end
end


[~,maxcovind] = max(cov');
%t = -0.2:0.01:0.2;
timebins = out.corr.time(out.bins);
maxcov= timebins(maxcovind);



t = -0.2:0.02:0.2;
%h = histcounts(maxcov,t,'Normalization','probability');
h = histcounts(maxcov,t);
figure
bar(t(2:end),h);

hold on
plot([0 0],[0 max(h)],'k--')

xlabel('Peak cross-covariance time (seconds)')
ylabel('Fraction of cell pairs');

figFile = [figDir,'Spiking\crosscov_CA1-PFC_selectivecells'];
print(figFile,'-djpeg');
print(figFile,'-dpdf');
end

%% --- Selective Cells Only - Same Odor Preference - Only look at cell pairs with same odor preference
if sameprefcells == 1
cov = [];
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
        
        eps = runeps(runeps(:,1) == day,2);
        
        for ep = eps'
            
            %find first good tet
            filt = '$numspikes > 0';
            tets = evaluatefilter(cellinfo{day}{ep},filt);
            tet = tets(1,1);
            
            eeg = loadeegstruct(animDir, animal, 'eeg',day,ep,tet);
            times = geteegtimes(eeg{day}{ep}{tet});
            
            
            wins = npWins{day}{ep};
            
            [wintime,winvec] = wb_list2vec(wins,times);
            nonnpvec = ~winvec;
            nonnplist = vec2list(nonnpvec,times);
            
            
            filt = 'strcmp($area,''CA1'') && strcmp($selectivity,''rightSelective'')';
            cellsCA1 = evaluatefilter(cellinfo{day}{ep},filt);
            
            
            filt = 'strcmp($area,''PFC'') && strcmp($selectivity,''leftSelective'')';
            cellsPFC = evaluatefilter(cellinfo{day}{ep},filt);
            
            [A,B] = meshgrid(1:size(cellsCA1,1),1:size(cellsPFC,1));
            pairInds=reshape(cat(2,A',B'),[],2);
            
            for p = 1:size(pairInds,1)
                
                ind = [day, ep, cellsCA1(pairInds(p,1),:), cellsPFC(pairInds(p,2),:)];
                
                out = DFAsj_getthetacrosscov(ind, nonnplist, spikes);
                
                if any(~isnan(out.Zcrosscov_sm)) && out.Neventscorr > 10
                cov = [cov; out.Zcrosscov_sm];
                end
            end
            
            filt = 'strcmp($area,''CA1'') && strcmp($selectivity,''leftSelective'')';
            cellsCA1 = evaluatefilter(cellinfo{day}{ep},filt);
            
            
            filt = 'strcmp($area,''PFC'') && strcmp($selectivity,''rightSelective'')';
            cellsPFC = evaluatefilter(cellinfo{day}{ep},filt);
            
            [A,B] = meshgrid(1:size(cellsCA1,1),1:size(cellsPFC,1));
            pairInds=reshape(cat(2,A',B'),[],2);
            
            for p = 1:size(pairInds,1)
                
                ind = [day, ep, cellsCA1(pairInds(p,1),:), cellsPFC(pairInds(p,2),:)];
                
                out = DFAsj_getthetacrosscov(ind, nonnplist, spikes);
                
                if any(~isnan(out.Zcrosscov_sm)) && out.Neventscorr > 10
                cov = [cov; out.Zcrosscov_sm];
                end
            end

            
        end
    end
end


[~,maxcovind] = max(cov');
%t = -0.2:0.01:0.2;
timebins = out.corr.time(out.bins);
maxcov= timebins(maxcovind);



t = -0.2:0.02:0.2;
%h = histcounts(maxcov,t,'Normalization','probability');
h = histcounts(maxcov,t);
figure
bar(t(2:end),h);

hold on
plot([0 0],[0 max(h)],'k--')

xlabel('Peak cross-covariance time (seconds)')
ylabel('Fraction of cell pairs');

figFile = [figDir,'Spiking\crosscov_CA1-PFC_selectivecells_oppositeodorpref'];
print(figFile,'-djpeg');
print(figFile,'-dpdf');
end
