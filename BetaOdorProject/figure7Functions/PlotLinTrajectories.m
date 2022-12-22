    %% PlotLinTrajectories

%{
The gist of this is that we want to get mean tuning curves for all the
oyutboud and inbound trajectories.  The way I have this now is that I have
segmented each 'run' as a run from one well to another. 
The stepwise operation for this is as follows:
1. choose the trajectories you want to use for a given run
2. get a good set of positional bins for each frame
3. for each contiguous block of time, get the positional occupancy of the
animal, and interp the appropriate spikes so you get nspikes per bin
4. I choose to smooth after merging spikes / bin by time/bin
5. smooth your dataset
6. use this as a template fo rthe decoder which will be a mean firing rate 
 spikes per sec, in each bin.  Then you can use those as the poisson means
 when you ask what the rate ought to be at a given time.  Generally
 speaking a cell will fire 1, 2, or 3 times in a given ripple or theta
 cycle, so your posterior rates ought to be between 0 and 3 ish.
%}



%% 1. get your trajectories
% choose your session
i=6;



lincoords=table2array(SuperRat(i).LinCoords);
coords2=SuperRat(i).tracking.data;
% lets start with the left trajectory
keepinds=lincoords(:,4)==3 & lincoords(:,5)==1;

linLeft=lincoords(keepinds,:);
% get all the time breaks (this is the last coord before the new set
breaks=find(diff(linLeft(:,1))>1);
epochs=[[linLeft(1); linLeft(breaks+1,1)] [linLeft(breaks); linLeft(end,1)]];
epochinds=[[1; breaks+1] [breaks; size(linLeft,1)]];

%% sanity check
figure;
sp(1)=subplot(1,3,1);
% first the trajectory
for ei=1:length(epochinds)
    plot(linLeft(epochinds(ei,1):epochinds(ei,2),2),linLeft(epochinds(ei,1):epochinds(ei,2),3));
    hold on;
end
set(gca,'ydir','reverse');
title('all left runs');
% now the speed
sp(2)=subplot(1,3,2);
for ei=1:length(breaks)-1
    plot(linLeft(epochinds(ei,1):epochinds(ei,2),1),linLeft(epochinds(ei,1):epochinds(ei,2),7));
    hold on;
end
title('all speeds for those runs')
xlabel('time'); ylabel('speed');
sp(3)=subplot(1,3,3);

% and the linearized position
for ei=1:length(breaks)-1
    plot(linLeft(epochinds(ei,1):epochinds(ei,2),1),linLeft(epochinds(ei,1):epochinds(ei,2),8));
    hold on;
end
xlabel('time'); ylabel('linearized position');
linkaxes(sp(2:3),'x');


%% okay lets measure some occupancies and neurons
% first lets do occupancy

% get number of ts in each bin
pixocc=accumarray(linLeft(:,8),1);
%[occupancy,bins]=histcounts(linLeft(:,5),1:100);
figure; plot(1:length(pixocc),log2(pixocc))
xlabel('linearized position (percentile of run)'); ylabel('time spent');

% now lets see if we can grab a unit
for un=1:length(SuperRat(i).units)
    % grab all the occupied pixels
    [occupancy,bins]=histcounts(linLeft(:,8),1:max(linLeft(:,8)));
    spikes=SuperRat(i).units(un).ts;
    % now epoch the spikes based on when the runs are
    Espikes=EpochCoords(spikes(:,1),epochs);
    % get a position for each spike
    spikepos=interp1(linLeft(:,1),linLeft(:,8),Espikes,'nearest');
    % get number of spikes per position
    [nspikes,bins]=histcounts(spikepos,bins);
    % now plot the two out
    figure; subplot(3,1,1);
    plot(bins(1:end-1),occupancy);
    subplot(3,1,2);
    plot(bins(1:end-1),nspikes);
    subplot(3,1,3);
    plot(bins(1:end-1),SmoothMat2(nspikes./occupancy,[5 0],2));
end


%% okay now we can use the above to make a pretty figure

% for a given rat we first plot out our trajectories
%savedir=uigetdir;

for i=1:length(SuperRat)
    if SuperRat(i).longTrack
        figure;
        subplot(2,2,1);
        trajinds=[3 1; 3 2; 1 3; 2 3];
        linpull=[4 5; 6 7; 4 5; 6 7];
        keepinds=[];
        colors=jet(20); colors=colors([1 5 15 20],:);
        allposplot=[];
        % this plots the trajectories BEFORE the slow parts are taken out
        
        for tr=1:4
            % first gather the four trajectories by nanning out all the other data
            % this is smoothed and % normalized for each run
            temptraj=table2array(SuperRat(i).LinCoords);
            
            % take account of the correct trajectories at the correct speed
            keepinds=temptraj(:,4)==trajinds(tr,1) & temptraj(:,5)==trajinds(tr,2);

            % cat the real data here (ts,x,y,linpos,lindist)
            allposplot=[allposplot; temptraj(keepinds,:)];
            % grab a temporary trajectory
            temptraj(~keepinds,:)=nan;
            % and plot it
            plot(temptraj(:,2),temptraj(:,3),'color',colors(tr,:));
            hold on;
        end
        box off; axis off;
        title([SuperRat(i).name ' Day # ' num2str(SuperRat(i).daynum)]);
        % and nan out the slow times (dx for linear pos)
        allposplot=sortrows(allposplot,1); % order chronologically
        %allposplot(tooslow,:)=[];
        
        % get dwell time
        [occupancy,bins]=histcounts(allposplot(:,8),1:max(allposplot(:,8)));
        occupancy=accumarray(allposplot(:,8),1,[max(allposplot(:,8)),1]);
        
        % find the breaks between each run to epoch our spikes
        breaks=find(diff(allposplot(:,1))>1);
        epochs=[[allposplot(1); allposplot(breaks+1,1)] [allposplot(breaks,1); allposplot(end,1)]];
        epochinds=[[1; breaks+1] [breaks; size(allposplot,1)]];
        
        % now grab our cells (just taky pyrs
        cellnotes={SuperRat(i).units.type};
        for cn=1:length(cellnotes)
            if isempty(cellnotes{cn})
                cellnotes(cn)={'bad'};
            end
        end
        
        unitcriterion=contains(cellnotes,'pyr');
        unitdata=SuperRat(i).units(unitcriterion);
        fullcurves=[]; curvepeaks=[];
        cellsort=[];
        for j=1:length(unitdata)
            spikes=unitdata(j).ts;
            % nowe we have linearize position plots to concatenate
            % now epoch the spikes based on when the runs are
            Espikes=EpochCoords(spikes(:,1),epochs);
            % snap each spike to its nearest position
            spikepos=interp1(allposplot(:,1),allposplot(:,8),Espikes,'nearest');
            % get number of spikes per position
            [nspikes]=accumarray(spikepos,1,[max(allposplot(:,8)),1]);
            %[nspikes,bins]=histcounts(spikepos,bins);
            % now save the grand total tuning curve
            fullcurves(j,:)=SmoothMat2(nspikes./(occupancy/30),[0 20],5);
            cellsort(j,1)=contains(unitdata(j).area,'PFC')+1;
            [curvepeaks(j),cellsort(j,2)]=max(fullcurves(j,:));
        end
        
        deadcells=nansum(fullcurves,2)<0.5;
        fullcurves(isnan(fullcurves) | isinf(fullcurves))=0;
        fullcurves(deadcells,:)=[]; 
        cellsort(deadcells,:)=[]; 
        curvepeaks(deadcells)=[];
        % gotta sort these guys out
        [sorted,sortind]=sortrows(cellsort);
        
        % and now i can sort the linearized plots from left, right, out, and back
        % now make a line plot for each of the four trajectories
        % or do a filled line plot
        
        %subplot(2,2,2); imagesc(zscore(fullcurves(sortind,:),1,2));
        regcolors=[1 0 0; 0 0 1];
        subplot(2,2,[2 4])
        clear sp
        for k=1:size(fullcurves,1)
            sp(k)=patch([1 1:size(fullcurves,2) size(fullcurves,2)], size(fullcurves,1)-k+0.9*[0 fullcurves(sortind(k),:)...
                /max(fullcurves(sortind(k),:)) 0],...
                regcolors(cellsort(sortind(k),1),:),'EdgeColor','none');
        end
        xlabel('Distance from center well'); ylabel('cell number');
        ylim([0 size(fullcurves,1)]);
        legend(sp([1 length(sp)]),{'HPC','PFC'}); title('All Directions & Routes');
        yyaxis right;
        set(gca,'YTick',0:length(curvepeaks(sortind)),'YTickLabel',(curvepeaks(sortind)));
        set(gca,'YLim',[-0.5 length(curvepeaks(sortind))-0.5]);
        % now in bottom left try the four perms of direction
        
        %
        smallpos=[9 10 13 14];
        fullcurves=[]; clear cellsort;
        for tr=1:4
            % first gather the four trajectories by nanning out all the other data
            temptraj=SuperRat(i).LinCoords;
            keepinds=temptraj(:,4)==trajinds(tr,1) & temptraj(:,5)==trajinds(tr,2);
            temptraj=temptraj(keepinds,:); % pull fav line, and rows
            
            % and nan out the slow times (dx for linear pos)
            temptraj=sortrows(temptraj,1); % order chronologically
            
            
            % find the epochs sou you can kill bad spikes
            breaks=find(diff(temptraj(:,1))>1);
            epochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1); temptraj(end,1)]];
            
            % get occupancy
            [occupancy,bins]=histcounts(temptraj(:,8),1:max(allposplot(:,8)));

            % now for each unit capture a mean tuning curve for each run
            for j=1:length(unitdata)
                spikes=unitdata(j).ts;
                Espikes=EpochCoords(spikes(:,1),epochs); % only pull spike ts
                spikepos=interp1(temptraj(:,1),temptraj(:,8),Espikes,'nearest'); % now interp to position
                [nspikes,bins]=histcounts(spikepos,bins); % get number of spikes per position
                % now save this tuning curve
                fullcurves(j,tr,:)=SmoothMat2(nspikes./occupancy,[5 0],2);
                % and also, add to that cell
                %/// not sure how i want to save out these data, i'll need a bunch
                %for calcuating ripples
                cellsort(j,1)=contains(unitdata(j).area,'PFC')+1;
                if tr==1, [~,cellsort(j,2)]=max(fullcurves(j,tr,:)); end
            end
        end
        
        deadcells=sum(sum(fullcurves(:,:,:),3),2)<0.5;
        fullcurves(deadcells,:,:)=[]; cellsort(deadcells,:)=[];
        fullcurves(isnan(fullcurves) | isinf(fullcurves))=0;
        
        smalltitles={'out left','out right','in left','in right'};
        for tr=1:4
            subplot(4,4,smallpos(tr));
            % if its the first trajectory, save out order, otherwise use first
            % trajectories order and plot
            if tr==1
                [sorted,sortorder]=sortrows(cellsort);
                imagesc(zscore(squeeze(fullcurves(sortorder,tr,:)),1,2));
                if size(sorted,1)>1
                    hold on; plot([0 100],[find(diff(sorted(:,1))~=0) find(diff(sorted(:,1))~=0)]+.5,'k','LineWidth',2)
                end
            else
                imagesc(zscore(squeeze(fullcurves(sortorder,tr,:)),1,2));
                
                if size(sorted,1)>1
                    hold on; plot([0 100],[find(diff(sorted(:,1))~=0) find(diff(sorted(:,1))~=0)]+.5,'k','LineWidth',2)
                end
                
            end
            title(smalltitles{tr});
        end
        
        %savefig(fullfile(savedir, [SuperRat(i).name ' Day # ' num2str(SuperRat(i).daynum)]));
    end
end
    
   
    
    
    
    
