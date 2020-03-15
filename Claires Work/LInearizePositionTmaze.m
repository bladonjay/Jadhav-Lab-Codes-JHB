%LinearizePositionTmaze
%{
this produces a 'lincoords field that has 10 columns and n video frames
rows


The rows are as follows:
1. the timestamp of the frame
2. the x coordinate
3. the y coordinate
4. the distance from the left linear trajectory
5. the index of the left trajectory
6. the distance from the right linear trajectory
7. the index of the right trajectory
8. the well the rat is leaving (1,2 outers, 3 home)
9. the well the rat will end up at
10. the epoch in the day that this is
11. the speed of the animal (I smoothed 7 because it oscillates HARD

The algorithm is as follows:
1. Outline the average trajectory for all unique trajectories
2. now for each position, grab the closest point on each trajectory
3. track the point, and the distance to that point

4. go back over all the trajectories, dice into blocks between when the
animal is in close proximity to a well
5. categorise the whole trajectory by which well it came from and which
well it went to
6. velocity filter and ripple filter




%}
% linearize position for claires data
%%
% first work off of 'superRat' struct
% at the bottom of this is where we designate the side well and home well
% distances

for i=1:length(SuperRat)
    tic
    coorddata=SuperRat(i).tracking.data;
    % these coord data come in as ts, x, y, dir?, smoothed vel (not
    % really), and epoch
    
    % first gather a typical route for these animals.
    useinds=ismember(coorddata(:,6),SuperRat(1).RunEpochs');
    if sum(useinds)>0
        
        if ~isfield(SuperRat(i),'mazeMap') || contains(SuperRat(i).name,'CS39')
            fprintf('No drawn tracks for this session, draw new ones: \n');
            fprintf('First do left track, then do right track \n');
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
            for tk=1:length(tracks)
                trackpos=[round(xpos(tracks==tk)*10) round(ypos(tracks==tk)*10)];
                % have to expand each within
                for j=1:size(trackpos,1)-1
                    % expand both and see which is larger
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
            
        else
            mytracks=SuperRat(i).mazeMap.mytracks;
        end
    end
    % now you need to identify runs, this will be based on the trials that are
    % in the data struct.  Basically, you ought to be able to go out from the
    % odor ends and project out till the next odor start
    
    % so to start, lets grab the first coord that puts the animal in the home
    % well.  This will be the first coord whose radius from the home well pixel
    % is within say 10% of the total x and y distance of the maze
    alltracks=[mytracks{1}; mytracks{2}];
    homewell.center=nanmean([mytracks{1}(1,:); mytracks{2}(1,:)],1);
    homewell.radius=max(hypot(alltracks(:,1),alltracks(:,2)))*.05; % the radius is designated here *********
    sidewells.center=[mytracks{1}(end,:); mytracks{2}(end,:)];
    sidewells.radius=homewell.radius;
    
    SuperRat(i).mazeMap.homewell=homewell;
    SuperRat(i).mazeMap.sidewells=sidewells;
    SuperRat(i).mazeMap.mytracks=mytracks;
    SuperRat(i).mazeMap.notes='Used trajectories to map these positions out, you may want to redo this';
    fprintf('Linearized Coords for Session %d %s day %d in %.2f minutes \n', i,SuperRat(i).name, SuperRat(i).daynum, toc/60);
end

%%

% now go until you have a radius small enough to start this is
% determinisitc
for i=15:length(SuperRat)
    if SuperRat(i).longTrack>0
        tic
        coorddata=SuperRat(i).tracking.data;
        AllLinCoords=[];
        homewell=SuperRat(i).mazeMap.homewell;
        sidewells=SuperRat(i).mazeMap.sidewells;
        mytracks=SuperRat(i).mazeMap.mytracks;
        for epoch=1:length(SuperRat(i).RunEpochs)
            % start with the behavior epochs
            thesecoords=coorddata(coorddata(:,6)==SuperRat(i).RunEpochs(epoch),:);
             
            % need to remove the oscillation
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
            % templincoords are: ts, x, y,    Ldist, Lind, Rdist,Rind, origin, destination
            tempLinCoords=[thesecoords(1:firststart-1,1:3) nan(firststart-1,8)];
            goaldists=nan(firststart-1,3); % this will be the distance of this timestamp to each goal
            
            
            
            % initialize the first datapoint here
            [Ldist,Lind]=min(hypot(abs(mytracks{1}(:,1)-thesecoords(thispos,2)),...
                abs(mytracks{1}(:,2)-thesecoords(thispos,3))));
            % left goal dist
            Lgdist=hypot(abs(sidewells.center(1,1)-thesecoords(thispos,2)),...
                abs(sidewells.center(1,2)-thesecoords(thispos,3)));
            
            % right route and goal dist and index
            [Rdist,Rind]=min(hypot(abs(mytracks{2}(:,1)-thesecoords(thispos,2)),...
                abs(mytracks{2}(:,2)-thesecoords(thispos,3))));
            Rgdist=hypot(abs(sidewells.center(2,1)-thesecoords(thispos,2)),...
                abs(sidewells.center(2,2)-thesecoords(thispos,3)));
            
            Hgdist=hypot(abs(homewell.center(1)-thesecoords(thispos,2)),...
                abs(homewell.center(2)-thesecoords(thispos,3)));
            
            
            % goal distances are the distances from each goal
            goaldists(thispos,:)=[Lgdist Rgdist Hgdist];
            % load up this coord val
            tempLinCoords(thispos,:)=[thesecoords(thispos,1:3) Ldist Lind Rdist Rind nan nan nan thisspeed(thispos)];
            thispos=thispos+1;
            
            % while we havent done all the coords...
            while thispos<size(thesecoords,1)
                
                % now while he hasnt reached an end arm yet
                while min(goaldists(end,:))>sidewells.radius % basically while hes not close to any other arm
                    
                    % now for this trajectory go until youre close to a reward well
                    % left min dist, and left linearized index of closest
                    [Ldist,Lind]=min(hypot(abs(mytracks{1}(:,1)-thesecoords(thispos,2)),...
                        abs(mytracks{1}(:,2)-thesecoords(thispos,3))));
                    % left goal dist
                    Lgdist=hypot(abs(sidewells.center(1,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(1,2)-thesecoords(thispos,3)));
                    
                    % right route and index and goal dist
                    [Rdist,Rind]=min(hypot(abs(mytracks{2}(:,1)-thesecoords(thispos,2)),...
                        abs(mytracks{2}(:,2)-thesecoords(thispos,3))));
                    Rgdist=hypot(abs(sidewells.center(2,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(2,2)-thesecoords(thispos,3)));
                    
                    % home goal distance
                    Hgdist=hypot(abs(homewell.center(1)-thesecoords(thispos,2)),...
                        abs(homewell.center(2)-thesecoords(thispos,3)));
                    % load up our coords
                    tempLinCoords(thispos,:)=[thesecoords(thispos,1:3) Ldist Lind Rdist Rind nan nan nan thisspeed(thispos)];
                    goaldists(thispos,:)=[Lgdist Rgdist Hgdist];
                    
                    % increment
                    thispos=thispos+1;
                    if thispos>size(thesecoords,1)
                        %keyboard % dbcont to make sure we're okay
                        break
                    end
                end
                
                % label this turn based on where he came from and where he's going,
                % generally even if the rat strays out to one side, the trajectories
                % will go from one well to another
                
                [~,destination]=min(goaldists(end,:)); % where is the end?
                % add our origin and destination to the matrix
                tempLinCoords(thisstart:thispos-1,8)=origin;
                tempLinCoords(thisstart:thispos-1,9)=destination;
                
                % and we'll save every lincoord so we can tell which is closest
                origin=destination; thisstart=thispos;
                
                % now wait till he leaves this well to retrack
                while min(goaldists(end,:))<sidewells.radius % basically while hes not close to any other arm
                    
                    % now for this trajectory go until youre close to a reward well
                    % left min dist, and left linearized index of closest
                    [Ldist,Lind]=min(hypot(abs(mytracks{1}(:,1)-thesecoords(thispos,2)),...
                        abs(mytracks{1}(:,2)-thesecoords(thispos,3))));
                    % left goal dist
                    Lgdist=hypot(abs(sidewells.center(1,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(1,2)-thesecoords(thispos,3)));
                    
                    % right route and index and goal dist
                    [Rdist,Rind]=min(hypot(abs(mytracks{2}(:,1)-thesecoords(thispos,2)),...
                        abs(mytracks{2}(:,2)-thesecoords(thispos,3))));
                    Rgdist=hypot(abs(sidewells.center(2,1)-thesecoords(thispos,2)),...
                        abs(sidewells.center(2,2)-thesecoords(thispos,3)));
                    
                    % home goal distance
                    Hgdist=hypot(abs(homewell.center(1)-thesecoords(thispos,2)),...
                        abs(homewell.center(2)-thesecoords(thispos,3)));
                    % load up our coords
                    tempLinCoords(thispos,:)=[thesecoords(thispos,1:3) Ldist Lind Rdist Rind nan nan nan thisspeed(thispos)];
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
                tempLinCoords(thisstart:thispos-1,8)=0; % to say he hasnt left anywhere
                tempLinCoords(thisstart:thispos-1,9)=destination; % destination is where hes at this time
                
                % and we'll save every lincoord so we can tell which is closest
                origin=destination; thisstart=thispos;
                
                
                % validate by showing that the correcty trajectory has a smaller sum of
                % differences,
                % and now assign that trajectory to all the previous points
            end
            AllLinCoords=[AllLinCoords; tempLinCoords];
        end
        fprintf('Session %d done in %.2f seconds \n',i,toc);
        SuperRat(i).LinCoords=AllLinCoords;
    else
        fprintf('Session %d was on a short track \n',i);
    end
end





%% now to designate the epoch for each video frame

for ses=1:length(SuperRat)
    myCoords=SuperRat(ses).LinCoords;
    allTracking=SuperRat(ses).tracking.data;
    % this will be the last datapoint until the next epoch starts
    epochbreaks=allTracking([find(diff(allTracking(:,6))); size(allTracking,1)],1);
    
    % something like start from 1:break(1), then from 2:end
    % (break(i)+1):break(i+1)
    myCoords(:,10)=1;
    for i=1:length(epochbreaks)-1
        myCoords(myCoords(:,1)>epochbreaks(i) & myCoords(:,1)<=epochbreaks(i+1),10)=i+1;
    end
    SuperRat(ses).LinCoords=myCoords;
end

%% now some sanity checks


figure;
subplot(6,1,1);















