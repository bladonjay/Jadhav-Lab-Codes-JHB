%% this is going to be a generalization of the Linearize position T maze script
%{
The idea is that you can generate any number of prototypical trajectories
to goal locations, then you can use this script to find each run from one
goal location to another, and plot the trajectory to that goal
% You can have multiple routes to each goal, and the algo will choose the
track that the animal is most often closest to.



%}
%%
% this is for two side wells, the first two tracks are for the outbound
% trajectories, and the next two tracks are for the inbound trajectories.


for i=1:length(SuperRat)
    tic
    if length(SuperRat(i).mazeMap.mytracks)<4
        coorddata=SuperRat(i).tracking.data;
        % these coord data come in as ts, x, y, dir?, smoothed vel (not
        % really), and epoch

        % first gather only the used epochs
        useinds=ismember(coorddata(:,6),SuperRat(1).RunEpochs');

        if sum(useinds)>0


            fprintf('Not enough drawn tracks for this session, draw new ones: \n');
            fprintf('First do left out, then rigth out, left back right back \n');
            [xpos,ypos,tracks]=outlinemaze2(coorddata(useinds,(1:3)));

            % now generate a reasonable number of continuous points that will grab
            % these trajectories

            % first round these data to a decent number
            % empricially, we can probably use really small bins, and then when we get
            % linpos we can expand it by smoothing across tiny bins, and then taking
            % larger bins.
            % first generate a grid of positions to snap to

            % now interpolate the tracks and scale up
            mytracks=cell(length(unique(tracks)),1);
            for tk=1:length(unique(tracks))
                % these are the coordinates
                trackpos=[round(xpos(tracks==tk)*10) round(ypos(tracks==tk)*10)];
                % have to expand each within
                for j=1:size(trackpos,1)-1
                    % expand both, basically interpolate by each step
                    % the step will always be a single pixel in one
                    % dimension
                    numbins=max(abs(diff(trackpos(j:j+1,:),1)));
                    trackx=round(linspace(trackpos(j,1),trackpos(j+1,1),numbins));
                    tracky=round(linspace(trackpos(j,2),trackpos(j+1,2),numbins));
                    mytracks{tk}=[mytracks{tk}; [trackx' tracky']];
                end
                % now convolve with the gaussian kernel that has a mean of say two bins
                mytracks{tk}(:,1)=SmoothMat2(mytracks{tk}(:,1),[0 25],5);
                mytracks{tk}(:,2)=SmoothMat2(mytracks{tk}(:,2),[0 25],5);

                % and finally bring back to normal coordinates
                mytracks{tk}=mytracks{tk}/10;
            end

        end
        % now you need to identify runs, this will be based on the trials that are
        % in the data struct.  Basically, you ought to be able to go out from the
        % odor ends and project out till the next odor start

        % so to start, lets grab the first coord that puts the animal in the home
        % well.  This will be the first coord whose radius from the home well pixel
        % is within say 10% of the total x and y distance of the maze
        try alltracks=cell2mat(mytracks); catch alltracks=cell2mat(mytracks'); end
        homewell.center=nanmean([mytracks{1}(1,:); mytracks{2}(1,:)],1); % take the mean of the two starts
        % take the point that is the farthest from the home well, and then go
        % to 5% of that distance
        homewell.radius=max(hypot(alltracks(:,1),alltracks(:,2)))*.05; % the radius is designated here *********
        sidewells.center=[mytracks{1}(end,:); mytracks{2}(end,:)];
        sidewells.radius=homewell.radius;

        SuperRat(i).mazeMap.homewell=homewell;
        SuperRat(i).mazeMap.sidewells=sidewells;
        SuperRat(i).mazeMap.mytracks=mytracks;
        SuperRat(i).mazeMap.notes='Used trajectories to map these positions out, you may want to redo this';
        fprintf('Linearized Coords for Session %d %s day %d in %.2f minutes \n', i,SuperRat(i).name, SuperRat(i).daynum, toc/60);

    end
end
%%
% now go until you have a radius small enough to start this is
% deterministic
% this creates fields AllLinCoords


%The cols are: 1. ts, 2. x, 3. y, 4. origin, 5. destination, 
% 6. epoch,  7. speed, 8. Lodist, 9. Loind, 10. Rodist, 11. Roind,
% 12 Lidist 13 Liind 14 ridist 15 riind



for i=1:length(SuperRat)
    if SuperRat(i).longTrack>0
        tic
        coorddata=SuperRat(i).tracking.data;
        AllLinCoords=[];
        homewell=SuperRat(i).mazeMap.homewell;
        sidewells=SuperRat(i).mazeMap.sidewells;
        mytracks=SuperRat(i).mazeMap.mytracks;
        if SuperRat(i).longTrack==0
            homewell.radius=homewell.radius*3;
            sidewells.radius=homewell.radius;
        end
        for epoch=1:length(SuperRat(i).RunEpochs)
            % start with the behavior epochs
            thesecoords=coorddata(coorddata(:,6)==SuperRat(i).RunEpochs(epoch),:);
             
            % need to remove the oscillation (due to the uneven sampling
            % rate of the camera)
            thisspeed=smoothdata(thesecoords(:,5));
            
            % grab first coord thats close to the home well this is where we start
            firststart=find(hypot(abs(thesecoords(:,2)-homewell.center(1)),...
                abs(thesecoords(:,3)-homewell.center(2)))<homewell.radius,1,'first');
            
            fprintf('LinCoords will start %.2f seconds from the start of this epoch \n',...
                thesecoords(firststart,1)-thesecoords(1));
            % okay start somewhere,
            thispos=firststart; thisstart=firststart;
            origin=3; % start with home well as origin
            
            % initialize a linearized coordinate pos and the distance on each
            % trajectory
            % templincoords are: ts, x, y, origin, destination, epoch, speed, Ldist, Lind, Rdist,Rind,
            tempLinCoords=[thesecoords(1:firststart-1,1:3) nan(firststart-1,4) nan(firststart-1,2*length(mytracks))];
            goaldists=nan(firststart-1,3); % this will be the distance of this timestamp to each goal
            
            
            
            % initialize the first datapoint here
            
            % left goal dist, right goal distance and home goal distance
            Lgdist=hypot(abs(sidewells.center(1,1)-thesecoords(thispos,2)),...
                abs(sidewells.center(1,2)-thesecoords(thispos,3)));
            
            Rgdist=hypot(abs(sidewells.center(2,1)-thesecoords(thispos,2)),...
                abs(sidewells.center(2,2)-thesecoords(thispos,3)));
            
            Hgdist=hypot(abs(homewell.center(1)-thesecoords(thispos,2)),...
                abs(homewell.center(2)-thesecoords(thispos,3)));
            
            % and the trajectory point distances
            tjdist=[]; tjind=[];
            for tj=1:length(mytracks)
                [tjdist(tj),tjind(tj)]=min(hypot(abs(mytracks{tj}(:,1)-thesecoords(thispos,2)),...
                    abs(mytracks{tj}(:,2)-thesecoords(thispos,3))));
            end
            goaldists(thispos,:)=[Lgdist Rgdist Hgdist];
            lintracks=linearize([tjind; tjdist]);
            % load up this coord val
            tempLinCoords(thispos,:)=[thesecoords(thispos,1:3) nan nan SuperRat(i).RunEpochs(epoch), ...
                thisspeed(thispos) lintracks'];
            thispos=thispos+1;
            
            % while we havent done all the coords...
            while thispos<size(thesecoords,1)
                
                % now while he hasnt reached an end arm yet
                while min(goaldists(end,:))>sidewells.radius % basically while hes not close to any other arm
                    
                    
                    % left goal dist
                    Lgdist=hypot(abs(sidewells.center(1,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(1,2)-thesecoords(thispos,3)));
                    
                    % right route and index and goal dist                 
                    Rgdist=hypot(abs(sidewells.center(2,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(2,2)-thesecoords(thispos,3)));
                    
                    % home goal distance
                    Hgdist=hypot(abs(homewell.center(1)-thesecoords(thispos,2)),...
                        abs(homewell.center(2)-thesecoords(thispos,3)));
                    
                    % and the trajectory point distances
                    tjdist=[]; tjind=[];
                    for tj=1:length(mytracks)
                        [tjdist(tj),tjind(tj)]=min(hypot(abs(mytracks{tj}(:,1)-thesecoords(thispos,2)),...
                            abs(mytracks{tj}(:,2)-thesecoords(thispos,3))));
                    end
                    goaldists(thispos,:)=[Lgdist Rgdist Hgdist];
                    lintracks=linearize([tjind; tjdist]);
                    % load up this coord val
                    tempLinCoords(thispos,:)=[thesecoords(thispos,1:3) nan nan SuperRat(i).RunEpochs(epoch), ...
                        thisspeed(thispos) lintracks'];

                    % load up our coords
                    goaldists(thispos,:)=[Lgdist Rgdist Hgdist];
                    
                    % increment
                    thispos=thispos+1;
                    if thispos>=size(thesecoords,1)
                        %keyboard % dbcont to make sure we're okay
                        break
                    end
                end
                
                % label this turn based on where he came from and where he's going,
                % generally even if the rat strays out to one side, the trajectories
                % will go from one well to another
                
                [~,destination]=min(goaldists(end,:)); % where is the end?
                % add our origin and destination to the matrix
                tempLinCoords(thisstart:thispos-1,4)=origin;
                tempLinCoords(thisstart:thispos-1,5)=destination;
                
                % and we'll save every lincoord so we can tell which is closest
                origin=destination; thisstart=thispos;
                
                % now wait till he leaves this well to retrack
                while min(goaldists(end,:))<sidewells.radius % basically while hes not close to any other arm
                    
                    % now for this trajectory go until youre close to a reward well
                    % left min dist, and left linearized index of closest
                   
                    % left goal dist
                    Lgdist=hypot(abs(sidewells.center(1,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(1,2)-thesecoords(thispos,3)));
                    
                    % right route and index and goal dist                  
                    Rgdist=hypot(abs(sidewells.center(2,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(2,2)-thesecoords(thispos,3)));
                    
                    % home goal distance
                    Hgdist=hypot(abs(homewell.center(1)-thesecoords(thispos,2)),...
                        abs(homewell.center(2)-thesecoords(thispos,3)));
                    % load up our coords
                    tjdist=[]; tjind=[];
                    for tj=1:length(mytracks)
                        [tjdist(tj),tjind(tj)]=min(hypot(abs(mytracks{tj}(:,1)-thesecoords(thispos,2)),...
                            abs(mytracks{tj}(:,2)-thesecoords(thispos,3))));
                    end
                    goaldists(thispos,:)=[Lgdist Rgdist Hgdist];
                    lintracks=linearize([tjind; tjdist]); % alternate pos, dist, pos dist for all the trajectories
                    % load up this coord val
                    tempLinCoords(thispos,:)=[thesecoords(thispos,1:3) nan nan SuperRat(i).RunEpochs(epoch), ...
                        thisspeed(thispos) lintracks'];
                    goaldists(thispos,:)=[Lgdist Rgdist Hgdist];
                    
                    % increment
                    thispos=thispos+1;
                    if thispos>size(thesecoords,1)
                        %keyboard % dbcont to make sure we're okay
                        break
                    end
                end
                [~,destination]=min(goaldists(end,:)); % where is the end?
                % add our origin and destination to the matrix
                tempLinCoords(thisstart:thispos-1,4)=0; % to say he hasnt left anywhere
                tempLinCoords(thisstart:thispos-1,5)=destination; % destination is where hes at this time
                
                % and we'll save every lincoord so we can tell which is closest
                origin=destination; thisstart=thispos;
                
                
                % validate by showing that the correcty trajectory has a smaller sum of
                % differences,
                % and now assign that trajectory to all the previous points
            end
            AllLinCoords=[AllLinCoords; tempLinCoords];
        end
        fprintf('Session %d done in %.2f seconds \n',i,toc);
        SuperRat(i).AllLinCoords=AllLinCoords;
        SuperRat(i).AllLinCoords=array2table(AllLinCoords,'VariableNames',...
            {'ts', 'x', 'y', 'origin','destination', 'epoch','speed','Lodist',...
            'Loind','Rodist', 'Roind','Lidist','Liind','Ridist','Riind'});
% 12 Lidist 13 Liind 14 ridist 15 riind'}
    else
        fprintf('Session %d was on a short track \n',i);
    end
end
%%
% and epoch each coord
for ses=1:length(SuperRat)
    if SuperRat(ses).longTrack==1
    myCoords=SuperRat(ses).AllLinCoords;
    allTracking=SuperRat(ses).tracking.data;
    % this will be the last datapoint until the next epoch starts
    epochbreaks=allTracking([find(diff(allTracking(:,6))); size(allTracking,1)],1);
    
    % something like start from 1:break(1), then from 2:end
    % (break(i)+1):break(i+1)
    myCoords.epoch(:)=1;
    for i=1:length(epochbreaks)-1
        myCoords.epoch(myCoords.ts>epochbreaks(i) & myCoords.ts<=epochbreaks(i+1))=i+1;
    end
    
    % and normalize all the trajectories
    NumTraj=(size(myCoords,2)-7)/2;
    trajNames=myCoords.Properties.VariableNames;
    for tr=9:2:7+NumTraj*2
        
        [~,~,tempcoords]=histcounts(myCoords.(trajNames{tr}),100); % bin into 100 increments
        tempcoords(tempcoords==0)=nan;
        myCoords.(trajNames{tr})=tempcoords;
    end
    
    SuperRat(ses).AllLinCoords=myCoords;

    else
        SuperRat(ses).AllLinCoords=[];
    end
end




%% this is where we distill all the trajectories into the most likely traj
verbose=0;
CoordVarNames={'ts','x','y','originWell','destinationWell','epoch','speed',...
    'prefTrajDist','prefTrajInd'};

for i=1:length(SuperRat)
    if SuperRat(i).longTrack==1
    allVarnames=SuperRat(i).AllLinCoords.Properties.VariableNames;

    % this will asssign each trajectory based on origin and destination
    trajinds=[3 1; 3 2; 1 3; 2 3];
    
    if size(SuperRat(i).AllLinCoords,2)>11
        linpull=[8 9; 10 11; 12 13; 14 15];
    else
        linpull=[8 9; 10 11; 8 9; 10 11];
    end
    keepinds=[];
    colors=jet(20); colors=colors([1 5 15 20],:);
    allposplot=[];
    % this plots the trajectories BEFORE the slow parts are taken out
    if verbose,    figure; end
    % also you have to normalize all the trajectories
    for tr=1:4
        % first gather the four trajectories by nanning out all the other data

        temptraj=table2array(SuperRat(i).AllLinCoords);
        % filter based on speed % we'll do 3... but we need
        temptraj(temptraj(:,7)<3,:)=[];
        
        % take account of the correct trajectories at the correct speed
        keepinds=temptraj(:,4)==trajinds(tr,1) & temptraj(:,5)==trajinds(tr,2);
        %tooslow=SmoothMat2(abs([0; diff(temptraj(:,linpull(tr,2)))]),[0 50],10)<1;
        % cat the real data here (ts,x,y,linpos,lindist)
        
        allposplot=[allposplot; temptraj(keepinds,[1:7 linpull(tr,:)])];
        % grab a temporary trajectory
        temptraj(~keepinds,:)=nan;
        trajplot{tr}=temptraj(:,[2 3]);
        if verbose
            % and plot it
            subplot(3,4,tr);
            plot(temptraj(:,2),temptraj(:,3),'color',colors(tr,:));
            subplot(3,4,tr+4);
            plot(temptraj(keepinds,linpull(tr,2)),temptraj(keepinds,linpull(tr,1)),'color',colors(tr,:));
            subplot(3,4,tr+8);
            [a,b]=histcounts(temptraj(keepinds,linpull(tr,1)),1:100);
            bar(a);
            title(sprintf('%s %d', SuperRat(i).name, SuperRat(i).daynum));
            set(gcf,'Position',[1000,270,560,1060])
        end
    end
    

    SuperRat(i).LinCoords=array2table(sortrows(allposplot,1),...
        'VariableNames',CoordVarNames);
    else
        SuperRat(i).LinCoords=[];
    end
end
%%
        
        
rowinds=1000:2500;
figure; subplot(5,1,1); scatter(AllLinCoords(rowinds,2),AllLinCoords(rowinds,3),15,rowinds,'filled');
sp=subplot(5,1,2);  plot(AllLinCoords(rowinds,2)); yyaxis right; plot(AllLinCoords(rowinds,3));
sp(2)=subplot(5,1,3); plot(AllLinCoords(rowinds,8)); yyaxis right; plot(AllLinCoords(rowinds,10));
sp(3)=subplot(5,1,4); plot(AllLinCoords(rowinds,12)); yyaxis right; plot(AllLinCoords(rowinds,14));
sp(4)=subplot(5,1,5); plot(AllLinCoords(rowinds,7)); hold on;
scatter(1:length(rowinds),ones(1,length(rowinds)),15,rowinds,'filled');
linkaxes(sp,'x')