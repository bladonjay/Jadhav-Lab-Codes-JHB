% LFP Spatial plots at bottom of figure 2


%% for beta power first
for i=1:length(SuperRat)

    % first get spatial bins as whole numbers for each position
    % so it looks like my last plot uses something like 80 bins, which is
    % like a 1 cm/bin.  so the animal goes from about 25 to 121 in the x,
    % i think its about 100 pixels across, lets turn that into 80 pixels
    % across
    if SuperRat(i).longTrack==1
        % the factor is 5/4 then, first discretize
        binedges=1:1.2:150;
        runOK=ismember(SuperRat(i).AllLinCoords.epoch,SuperRat(i).RunEpochs);
        newx=discretize(SuperRat(i).AllLinCoords.x(runOK),binedges);
        newy=discretize(SuperRat(i).AllLinCoords.y(runOK),binedges);
        badcoords=isnan(newx) | isnan(newy);
        newx(badcoords)=[]; newy(badcoords)=[];
        % should be able to accumarray at mean all amplitudes for eacch spatial
        % bin- insert a nx2 matrix for first

        occup=accumarray([newx newy],1,[max(newx),max(newy)]);
        % probably smooth after dividing?
        figure; subplot(1,3,1);
        imagesc(SmoothMat2(occup,[10 10],2));
        set(gca,'CLim',[0 75]);
        title('occupancy')

        % need to downsample amp data to camera data, probably accumarray at
        % average getting average amp in high freq times for each low freq
        % interval
        slowtime=discretize(SuperRat(i).PFCbeta(:,1),SuperRat(i).AllLinCoords.ts(runOK));
        % for filtered, 1 ts, 2 filtered, 3 phase, 4 envelope
        usedata=~isnan(slowtime);
        
        % this is the mean amplitude over each 30 hz interval
        PFCampdata=accumarray(slowtime(usedata),SuperRat(i).PFCbeta(usedata,4),[],@mean);
        CA1ampdata=accumarray(slowtime(usedata),SuperRat(i).CA1beta(usedata,4),[],@mean);
        PFCampdata(badcoords)=[]; CA1ampdata(badcoords)=[];
        % this is those mean amplitudes put into spatial bins and averaged
        % again
        PFCamps=accumarray([newx(2:end) newy(2:end)],PFCampdata,[max(newx),max(newy)],@mean);
        CA1amps=accumarray([newx(2:end) newy(2:end)],CA1ampdata,[max(newx),max(newy)],@mean);

        % need to clean this up a tad... or just average over lots of sessions
        subplot(1,3,2);
        imagesc(SmoothMat2(PFCamps,[10 10],2));
        set(gca,'CLim',[0 75]);
        title('Mean PFC Beta Power');

        subplot(1,3,3);
        imagesc(SmoothMat2(CA1amps,[10 10],2));
        set(gca,'CLim',[0 75]);
        title('Mean CA1 Beta Power');
        sgtitle(sprintf('%s session %d', SuperRat(i).name, SuperRat(i).daynum))

        % both unsmoothed, you can fix later
        SuperRat(i).mazeMap.occup=occup;
        SuperRat(i).mazeMap.meanPFCBeta=PFCamps;
        SuperRat(i).mazeMap.meanCA1Beta=CA1amps;

        % maybe smooth, and then set any nonoccupied pixel back to zero
    end
end
    
    % n kill amps in which occupancy is too low

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