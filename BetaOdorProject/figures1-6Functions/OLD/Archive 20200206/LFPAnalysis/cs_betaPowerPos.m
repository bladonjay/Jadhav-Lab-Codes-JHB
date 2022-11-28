%plot position
%plot z-scored power

%cs_betaPowerPos('CS33', 3,2) 
function cs_betaPowerPos(animal, day,ep)
close all
[topDir,figDir] = cs_setPaths();
animDir = [topDir, animal,'Expt\',animal,'_direct\'];

regions = {'CA1','PFC','OB'};

pos = loaddatastruct(animDir, animal,'pos',day);
postime = pos{day}{ep}.data(:,1);
posdata = pos{day}{ep}.data(:,[2 3]);
for r = 1:length(regions)
    region = regions{r};
    
    tet = cs_getMostCellsTet(animal,day,ep,region);
    beta = loadeegstruct(animDir, animal, 'beta',day,ep,tet);
    times = geteegtimes(beta{day}{ep}{tet});
    timeinds = find(times>=postime(1) & times<=postime(end));
    data = double(beta{day}{ep}{tet}.data(:,3));
    data = data(timeinds);
    times = times(timeinds);
    
    %align pos time to eeg time
    newdata = interp1(times, data, postime);
    newdata(isnan(newdata)) = 0;
    
    %outbound traj only
        linstate = DFTFsj_getlinstate(animDir,animal,[day, ep], 6);        
        excludeState = getExcludePeriods(postime,(linstate{day}{ep}.state == 1 | linstate{day}{ep}.state == 3));
        goodtimeinds = ~isExcluded(postime, excludeState);
        
        newdata = newdata(goodtimeinds);
        posdata = posdata(goodtimeinds,:);

        
    %remove outliers before zscoring
    %outlierinds = find(isoutlier(newdata));
%     newdata(outlierinds) = [];
%     
%     %zscore
    newdata = (newdata- mean(newdata))/std(newdata);
%     
%     postime(outlierinds) = [];
%     posdata(outlierinds,:) = [];
    
    newpos = round(posdata);
    meanPower = accumarray(newpos,newdata',[],@mean);
    %imagesc(meanPower)
    
    meanPower(isnan(meanPower)) = 0;
    g = gaussian2(3,(3*3));
    smoothed = filter2(g, meanPower);
    imagesc(smoothed)
    
    [occupancy] = hist2(newpos(:,1), newpos(:,2), [1:size(meanPower,1)], [1:size(meanPower,2)]);
    
    smoothedocc = filter2(g, occupancy');
    lowocc = smoothedocc <= 1;
    
    smoothed(lowocc) = NaN;
    %meanPower = interp2(meanPower,3);
    %figure
    
    imagesc(smoothed');
    colormap(hot)
    title(region)
    colorbar
    
    figfile = [figDir,'Specgrams\',region,'_betaPower_position'];
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
end

