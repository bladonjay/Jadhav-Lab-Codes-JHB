%use CS33 day 3 ep 2
 
function cs_betaCoherencePos(animal, day,ep)
close all
[topDir, figDir] = cs_setPaths();
animDir = [topDir, animal,'Expt\',animal,'_direct\'];

regions = {'CA1-PFC','CA1-OB','PFC-OB'};

g = gaussian2(3,(3*3)); %for smoothing
bandpass = [15 30];


daystr = getTwoDigitNumber(day);

pos = loaddatastruct(animDir, animal,'pos',day);
postime = pos{day}{ep}.data(:,1);
posdata = pos{day}{ep}.data(:,[2 3]);

for r = 1:length(regions)
    region = regions{r};
    
    load([animDir,animal,'coherence',region,daystr]);
    times = coherence{day}{ep}.time;
    data = coherence{day}{ep}.Coh;
    goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
    data = data(find(goodrows),:);
    data = mean(data,1);
     
    timeinds = find(times>=postime(1) & times<=postime(end));
    data = data(timeinds);
    times = times(timeinds);
    
    %align pos time to eeg time
    newpos = interp1(postime, posdata, times);
%newpos(isnan(newpos)) = 0;
    %remove outliers before zscoring
    %outlierinds = find(isoutlier(data));
    %data(outlierinds) = [];
    
    
    %times(outlierinds) = [];
    %newpos(outlierinds,:) = [];
    newpos = round(newpos);
    meanPower = accumarray(newpos,data',[],@mean);
    smoothed = filter2(g, meanPower);
    
    [occupancy] = hist2(newpos(:,1), newpos(:,2), [1:size(meanPower,1)], [1:size(meanPower,2)]);
    
    smoothedocc = filter2(g, occupancy');
    lowocc = smoothedocc <= 1;
    
    smoothed(lowocc) = NaN;
    
%     figure
%     imagesc(meanPower);
    
    
    
    %remove bins where occupancy was low
    smoothed(lowocc) = -1;
    goodocc = smoothed(~lowocc);
    low = min(goodocc);
    high = max(goodocc);
    
    figure
    imagesc(smoothed');
    colormap(jet)
    title(region)
    colorbar
    figfile = [figDir,'Cohgrams\',region,'_betaCoherence_position'];
    
    saveas(gcf,figfile);
    print('-djpeg', figfile);
    print('-dpdf', figfile);
end
end
