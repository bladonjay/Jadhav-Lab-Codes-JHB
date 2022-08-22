% LFP Spatial plots at bottom of figure 2


%% for beta power first
for i=15:length(SuperRat)

    % first get spatial bins as whole numbers for each position
    % so it looks like my last plot uses something like 80 bins, which is
    % like a 1 cm/bin.  so the animal goes from about 25 to 121 in the x,
    % i think its about 100 pixels across, lets turn that into 80 pixels
    % across
    if ~isempty(SuperRat(i).PFCeegFile) && ~isempty(SuperRat(i).CA1eegFile) && SuperRat(i).longTrack==1
        % the factor is 5/4 then, first discretize
        binedges=1:.75:150;
        % using lincoords
        %{
        runOK=ismember(SuperRat(i).AllLinCoords.epoch,SuperRat(i).RunEpochs);
        newx=discretize(SuperRat(i).AllLinCoords.x(runOK),binedges);
        newy=discretize(SuperRat(i).AllLinCoords.y(runOK),binedges);
        %}
        % using all tracking
        runOK=ismember(SuperRat(i).tracking.data(:,6),SuperRat(i).RunEpochs);
        newx=discretize(SuperRat(i).tracking.data(runOK,2),binedges);
        newy=discretize(SuperRat(i).tracking.data(runOK,3),binedges);

        badcoords=isnan(newx) | isnan(newy);
        newx(badcoords)=[]; newy(badcoords)=[];
        % should be able to accumarray at mean all amplitudes for eacch spatial
        % bin- insert a nx2 matrix for first

        occup=accumarray([newx newy],1,[max(newx),max(newy)]);
        % probably smooth after dividing?
        figure; subplot(1,4,1);
        imagesc(SmoothMat2(occup,[10 10],2)');
        set(gca,'CLim',[0 75],'ydir','normal');
        title('occupancy')

        % need to downsample amp data to camera data, probably accumarray at
        % average getting average amp in high freq times for each low freq
        % interval
        % now load lfp
        
        load(SuperRat(i).PFCeegFile);
        
        % slowtime are the indices that match the timestmaps in tracking
        slowtime=discretize(betacontinuous(:,1),SuperRat(i).tracking.data(runOK,1));
        % for filtered, 1 ts, 2 filtered, 3 phase, 4 envelope
        usedata=~isnan(slowtime);
        
        % slowtime is the slow time indices of the fast data

        % this is the mean amplitude over each 30 hz interval
        PFCampdata=accumarray(slowtime(usedata),betacontinuous(usedata,4),[],@mean);
        PFCraw=betacontinuous(usedata,2);
        
        load(SuperRat(i).CA1eegFile);

        CA1ampdata=accumarray(slowtime(usedata),betacontinuous(usedata,4),[],@mean);
        CA1raw=betacontinuous(usedata,2);
        % the time data are the trackingdata (runok)
        
        PFCampdata(badcoords)=[]; CA1ampdata(badcoords)=[];
        % this is those mean amplitudes put into spatial bins and averaged
        % again
        PFCamps=accumarray([newx(2:end) newy(2:end)],PFCampdata,...
            [length(binedges),length(binedges)],@mean);
        CA1amps=accumarray([newx(2:end) newy(2:end)],CA1ampdata,...
            [length(binedges),length(binedges)],@mean);

        % need to clean this up a tad... or just average over lots of sessions
        subplot(1,4,2);
        imagesc(SmoothMat2(PFCamps,[10 10],2)');
        set(gca,'CLim',[0 75],'YDir','Normal'); ylim([0 140])
        title('Mean PFC Beta Power');

        subplot(1,4,3);
        imagesc(SmoothMat2(CA1amps,[10 10],2)');
        set(gca,'CLim',[0 75],'YDir','Normal'); ylim([0 140])
        title('Mean CA1 Beta Power');

        
        params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[20 30],'trialave',1);
        % and now coherency
        % window size 100 msec, step of 10 msec (150 and 30 samples)
        [C,~,~,~,~,t,f]=cohgramc(PFCraw,CA1raw,[1 .02],params);
        % turn our data in to real time by using our high sampling slow
        
        % already have the correct indices in slowtime, so now i have to
        % downsample those data into the length of C
        Cinds=linspace(0,length(C),length(slowtime(usedata))); Cinds(1)=1;
        % and then turn those inds into real data
        Cohinds=accumarray(ceil(Cinds)',slowtime(usedata)',[],@mean);

        % will make slightly more bins, but it'll be super close
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % now we have to snap the coherence data into camera time
        Cohampdata=accumarray(ceil(Cohinds), mean(C,2), [],  @mean);
        Cohampdata(badcoords)=[];
        % now we have Ctime and Cohampdata so we can snap it to space
        Cohamps=accumarray([newx(2:end) newy(2:end)],...
            Cohampdata,[length(binedges),length(binedges)],@mean);

        subplot(1,4,4);
        imagesc(SmoothMat2(Cohamps,[10 10],2)');
        set(gca,'CLim',[0 .8],'YDir','Normal'); ylim([0 140]);
        title('Mean CA1 Beta Power');
        sgtitle(sprintf('%s session %d', SuperRat(i).name, SuperRat(i).daynum))
        % both unsmoothed, you can fix later
        SuperRat(i).mazeMap.occup=occup;
        SuperRat(i).mazeMap.meanPFCBeta=PFCamps;
        SuperRat(i).mazeMap.meanCA1Beta=CA1amps;
        SuperRat(i).mazeMap.meanCA1PFCBeta=Cohamps;
        
        colormap([.4 .4 .4; hot(255)]);
        drawnow;

        % maybe smooth, and then set any nonoccupied pixel back to zero
    end
end
fprintf('\n got through it, saving now \n')
save('E:\Brandeis datasets\Claire Data\ClaireData2022-07-18','SuperRat','-v7.3')
    
    % n kill amps in which occupancy is too low
%%

% this is just testing the coherence data out
timeinds=SuperRat(i).tracking.data(runOK,1);
figure; subplot(3,1,1); plot(timeinds(1:20000),PFCampdata(1:20000));
hold on;
plot(SuperRat(i).trialdata.sniffstart,zeros(1,length(SuperRat(i).trialdata.sniffstart)),...
    'r.','MarkerSize',23)
subplot(3,1,2); plot(timeinds(1:20000),CA1ampdata(1:20000));
subplot(3,1,3); plot(timeinds(1:20000),Cohampdata(1:20000));
linkaxes(get(gcf,'Children'),'x')
%{

 %this is code from place plot to get a gist of the algorithm
% now snap coords to pixels to get occupancy
    
    noccup=accumarray([xpixels ypixels],velcoords(:,5),[max(xpixels),max(ypixels)]);
    % this sums up all the dt times for each coordinate
    % get noccup to be real time;
    
    
    
    % now snap spikes to their appropriate pixel
    
    velcoords(:,2)=xpixels; velcoords(:,3)=ypixels;
    
    newedges={1:max(xpixels),1:max(ypixels)};
    
    
    % now rebuild our cellspikes for Coords
    cellspikesC = interp1(velcoords(:,1),velcoords(:,1:end),newSpikeTimes,'nearest');
    % same filter
    cellspikesC = cellspikesC(cellspikes(:,4) > minvelocity,:);
    
% getting the number of spikes in each bin
[nspikes,~]=hist3(cellspikesC(:,2:3),'Edges',newedges);

%}

%% now plot out animal means and across animals
longsessions=[SuperRat.longTrack];
ratIDs=unique({SuperRat(longsessions==1).name});
mymap=[.4 .4 .4; hot(255)];

for i=1:length(ratIDs)
    fullPFCMap=[]; fullCA1Map=[]; fullCohMap=[];
    useSessions=find(contains({SuperRat.name},ratIDs{i}));
    for j=1:length(useSessions)
        % zscore to standardize across sessions
        fullPFCMap(:,:,j)=zscore(SuperRat(useSessions(j)).mazeMap.meanPFCBeta(:,1:145)',[],'all');
        fullCA1Map(:,:,j)=zscore(SuperRat(useSessions(j)).mazeMap.meanCA1Beta(:,1:145)',[],'all');
        fullCohMap(:,:,j)=zscore(SuperRat(useSessions(j)).mazeMap.meanCA1PFCBeta(:,1:145)',[],'all');
    end
    %figure;
    subplot(5,3,(i-1)*3+1);
    imagesc(SmoothMat2(mean(fullPFCMap,3),[10 10],1)); colormap(gca,mymap);
    title('PFC Beta Power'); set(gca,'YDir','normal'); axis off; box off;
    subplot(5,3,(i-1)*3+2);
    imagesc(SmoothMat2(mean(fullCA1Map,3),[10 10],1)); colormap(gca,mymap);
    title('CA1 Beta Power'); set(gca,'YDir','normal'); axis off; box off;
    subplot(5,3,(i-1)*3+3);
    imagesc(SmoothMat2(mean(fullCohMap,3),[10 10],1)); colormap(gca,mymap(:,[2 3 1]));
    title('CA1-PFC Beta Coherence');  set(gca,'YDir','normal'); axis off; box off;
    sgtitle(ratIDs{i})
end

 %%
% and now grand mean
fullPFCMap=[]; fullCA1Map=[]; fullCA1PFCMap=[];
longsessions=find([SuperRat.longTrack]);
mymap=[.4 .4 .4; hot(255)];
for i=1:length(longsessions)
    % need to zscore values after omitting zeros because its not working...
    oldmap=SuperRat(longsessions(i)).mazeMap.meanPFCBeta;
    oldmean=mean(linearize(oldmap(oldmap>0))); oldstd=std(linearize(oldmap(oldmap>0)));
    newmap=(oldmap-oldmean)./oldstd; newmap(oldmap==0)=nan;
    fullPFCMap(:,:,i)=newmap;
    oldmap=SuperRat(longsessions(i)).mazeMap.meanCA1Beta;
    oldmean=mean(linearize(oldmap(oldmap>0))); oldstd=std(linearize(oldmap(oldmap>0)));
    newmap=(oldmap-oldmean)/oldstd; newmap(oldmap==0)=nan;
    fullCA1Map(:,:,i)=newmap;
    oldmap=SuperRat(longsessions(i)).mazeMap.meanCA1PFCBeta;
    oldmean=mean(linearize(oldmap(oldmap>0))); oldstd=std(linearize(oldmap(oldmap>0)));
    newmap=(oldmap-oldmean)/oldstd; newmap(oldmap==0)=nan;
    fullCA1PFCMap(:,:,i)=newmap;
end
figure;
subplot(1,3,1);
imagesc(SmoothMat2(mean(fullPFCMap,3,'omitnan'),[10 10],2)');
title('PFC Beta Power');  set(gca,'CLim', [-.5 3], 'ydir','normal','Ylim',[0 150])
subplot(1,3,2);
imagesc(SmoothMat2(mean(fullCA1Map,3,'omitnan'),[10 10],2)'); colormap(mymap);
title('CA1 Beta Power');   set(gca,'CLim', [-.5 3], 'ydir','normal','Ylim',[0 150])
subplot(1,3,3);
imagesc(SmoothMat2(mean(fullCA1PFCMap,3,'omitnan'),[10 10],2)'); colormap(gca,mymap(:,[2 3 1]));
title('CA1-PFC Beta Coherence');  set(gca,'CLim', [-.5 3], 'ydir','normal','Ylim',[0 150])
sgtitle('All Rats')
%% then beta coherence
% will have ot calculate coherence for each session first

%% then ripples

% will have to identify ripples first
%{
ripple id algorithm per shantanu:
1. smooth amplitude using a gaussian with 4 msec std, go out wide
2.  min duration is 15 msec (min suprathresh), get that in samples
3. zscore the smoothed data
4. binarize for above ~2 or 4 sd above mean
5. extract contiguous bits where thresh is crossed
6. get starttime, endtime, mean amplitude, peak amplitude
%}