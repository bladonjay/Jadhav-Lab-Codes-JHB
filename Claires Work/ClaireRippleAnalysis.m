% first code to crossref the object coding with the place fields






%% okay lets show the mean maze coverage of the odor selective cells

% need more cells ot do this math...
ses=12;
odorcells=find(cellfun(@(a) ~isnan(a(1)) & a(2)<.05, {SuperRat(ses).units.OdorSelective}));
odorpref=cellfun(@(a) a(1)>a(2), {SuperRat(ses).units.OdorSelective});


lefttraj=zeros(2,100); leftcells=1;
righttraj=zeros(2,100); rightcells=1;
% now lets do some decoding
for i=1:length(odorcells)
    if ~isempty(SuperRat(ses).units(odorcells(i)).LinPlaceFields)
        if odorpref(odorcells(i))
            lefttraj=lefttraj+zscore(SuperRat(ses).units(odorcells(i)).LinPlaceFields);
            leftcells=leftcells+1;
        else
            righttraj=righttraj+zscore(SuperRat(ses).units(odorcells(i)).LinPlaceFields);
            rightcells=rightcells+1;
        end
    end
end
lefttraj=lefttraj/leftcells; righttraj=righttraj/rightcells;

figure; sp(1)=subplot(1,2,1); plot(lefttraj');
sp(2)=subplot(1,2,2); plot(righttraj'); linkaxes;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% RIPPLE DECODING NOW %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%
%{
General outline:
1. pick session
2. pick a ripple
3. make sure there are spikes in that ripple
3. parse ripple into 100ms timebins
4. for each hippocampal pyram, calculate the probability of each object, and of each
position
5. multiply products for each position together, and for each object
6. sum across positions is 1, sum across odors is 1
7. plot this out for each pos and each obj at each time
%}

% first lets get our ripples, we'll get the ripple, and then get which
% tetrodes it was on, and then get its duration

%
%% aggregate the ripples
%{
 so the ripples are discoordinated... super discoordinated.
I think i want to take the longest detected window, and then count the
number of tetrodes that detected any time within that window.
    
            
%}

ses=16;
% for now just going to choose whichever tetrode has the most valid units
% m units by n timestamps

% spatial prior (left) first get all the cells
PlacePriorL=cell2mat(cellfun(@(a) a(1,:), {SuperRat(ses).units.LinPlaceFields},'UniformOutput',false)');
PlacePriorR=cell2mat(cellfun(@(a) a(2,:), {SuperRat(ses).units.LinPlaceFields},'UniformOutput',false)');

% now just our units
HPCunits=contains({SuperRat(ses).units.area},'CA1') & contains({SuperRat(ses).units.type},'pyr') & ...
    (sum(PlacePriorL,2)>0 | sum(PlacePriorR,2)>0);
% object prior (right and left) need to get odor ids tho
ObjPrior=cell2mat({SuperRat(ses).units.OdorRates});

% same, threshold the units



[~,winningtet]=max(accumarray([SuperRat(ses).units.tet]',1));
winningrips=SuperRat(ses).ripdata([SuperRat(ses).ripdata.tetrode]==winningtet);

% make sure you have only one winning rip dataset
riptimes=[winningrips.starttime winningrips.endtime];
riptimes(diff(riptimes,1,2)<0.05,:)=[];

% okay lets grab a late ripple, say the 1500th one
for ripnum=350:380
    
    riptime=riptimes(ripnum,:);
    
    % validate that there are lotsa spikes
    HPCinds=find(HPCunits);
    clear membership;
    for j=1:sum(HPCunits)
        membership(j)=any(SuperRat(ses).units(HPCinds(j)).ts(:,1)>riptime(1) & SuperRat(ses).units(HPCinds(j)).ts(:,1)<riptime(2));
    end
    sum(membership)
    % okay lets dice this up into segments
    if sum(membership)>=5
        
        ripbins=riptime(1):.01:riptime(2);
        % in real nspikes (has to be poisson)
        % this is m units by n ripple bins, this is /100 cause that makes
        % its n spikes
        ripplemat=eventspikematrix(ripbins',SuperRat(ses).units(HPCunits),0,.01)'/100;
        
        % now to decode
        figure;
        subplot(2,2,1); imagesc([PlacePriorL; PlacePriorR]); hold on;
        plot([0 100],[length(PlacePriorL) length(PlacePriorL)],'k','LineWidth',3);
        subplot(2,2,2); imagesc(ripplemat(~dontuse,:));
        
        decodedL=[]; decodedR=[];
        for i=1:length(ripbins)
            % this is the vector for this ripple bin
            ripslice=ripplemat(~dontuse,i);
            
            
            % probmat will be a m trial, x n unit x t timebin
            
            % s is spike mat, fi=rate this unit at pos, ni=spikes this unit
            % p(pos|s)=c* prod across units i=1:n  (  p(pos)* fi(pos)^ni / fact(ni)  * exp(fi(pos)))
            % calculate the probability of each time given this activity
            % priorsmat is a m units by n spatial bins
            % that gets raised to the nspikes power fo reach position
            
            % second line is
            
            %probmat=(priorsmat.^repmat(ripslice,1,size(priorsmat,2)))./...
            %       repmat(factorial(ripslice),1,size(priorsmat,2)).*exp(-priorsmat);
            
            % tau is 0.01 (because nspikes is per 100 msec)
            
            %(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
            % prodmat= prod units (mean rate_i).^spikes_i)
            % summat= sum(mean_rate_i)
            
            probmatL=(PlacePriorL.^repmat(ripslice,1,size(PlacePriorL,2))).*...
                exp(-0.01*repmat(nanmean(PlacePriorL),size(PlacePriorL,1),1));
            % sum each units prob to 1
            probmatL=probmatL./sum(probmatL,2); % across time
            decodedL(i,:)=prod(probmatL);
            
            % Right side
            probmatR=(PlacePriorR.^repmat(ripslice,1,size(PlacePriorR,2))).*...
                exp(-0.01*repmat(nanmean(PlacePriorR),size(PlacePriorR,1),1));
            probmatR=probmatL./sum(probmatR,2); % across time
            decodedR(i,:)=prod(probmatR);
            
            
            % odor preference
            probmatL=(ObjPrior.^repmat(ripslice,1,size(ObjPrior,2))).*...
                exp(-0.01*repmat(nanmean(ObjPrior),size(ObjPrior,1),1));
            
        end
        
        decoded2=decodedL./nansum(decodedL); % summ all to 1
        subplot(2,2,3); imagesc(log2(decoded2)');
        set(gca,'CLim',[-5 0])
        %set(gca,'CLim',[-10 -6]);
        
        % now plot the object probabilities here
        
    end
    
    
end

%% do odor selective cells that code the same odor reactivate together?

% so you can get which ripples each cell is a member of and then
% crossreference that to odor coding.

% the problem though is that cells dont really code for odors that well

% first lets grab the vector of ripple starts and ends


% I think we'll want to plot the absolute sum of the two.  2 is
% perfectly coherent selectivity, 1 is two with exacyly half
% selectivity, or one who is positive strong selective, other weak
% negative, below one is two with opposing strong selectivities


it=1;
odorpairs=[];
for ses=1:length(SuperRat)
    % okay she doesnt have all the ripple values for all the tetrodes...
    % so sweep starting w most units until we get rips
    winningtet=accumarray([SuperRat(ses).units.tet]',1);
    [~,tetinds]=sort(winningtet,'descend');
    ripind=1;
    winningrips=[];
    fprintf('finding a suitable ripple wire \n');
    while isempty(winningrips)
        winningrips=SuperRat(ses).ripdata([SuperRat(ses).ripdata.tetrode]==tetinds(ripind));
        ripind=ripind+1;
    end
    fprintf('found good rip wire, running units \n');
    if ~isfield(SuperRat(ses),'ripcorrs')
        % make sure you have only one winning rip dataset
        riptimes=[winningrips.starttime winningrips.endtime];
        riptimes(diff(riptimes,1,2)<0.05,:)=[];
        
        %*** probably need to accompany the ripmat with the epoch each rip occured
        ripmat=eventspikematrix(riptimes(:,1),SuperRat(ses),0, riptimes(:,2)-riptimes(:,1));
        % only take ripples with 5 or more units participating
        ripmat(sum(ripmat>0,2)<=4,:)=[];
        % now see which code objects
        ripcorrs=corr(ripmat);
        SuperRat(ses).ripcorrs=ripcorrs;
    else
        ripcorrs=SuperRat(ses).ripcorrs;
    end
    % now gather the odor selectivity
    odorpref=[];
    for i=1:length(SuperRat(ses).units)
        % is cell odor selective?
        thispref=SuperRat(ses).units(i).OdorSelective;
        Rpref=diff(SuperRat(ses).units(i).OdorRates,[],2)>0; % reverse preference if needed
        thispref(Rpref,1)=-thispref(Rpref,1);
        % 4 is the cell number
        thispref(:,4)=i;
        if ~contains(SuperRat(ses).units(i).type,'pyr')
            thispref=nan(size(thispref,1),size(thispref,2));
        end
        odorpref=[odorpref; thispref];
    end
    %{
    % now gather odor preference
    odorpref=[];
    for i=1:length(SuperRat(ses).units)
        % pick best preference for this unit
        thispref=SuperRat(ses).units(i).OdorSelective; thispref(:,4)=i;
        Rpref=diff(SuperRat(ses).units(i).OdorRates(:,[1 2]),[],2)>0; % reverse preference if needed
        % replace selectivity score with which odor it preferred
        thispref= [Rpref thispref(:,2:4)];
        [a,bestind]=min(thispref(:,2)); % find smallest p
        if isnan(a), bestind=1; end % if there is no good p, use first ind (doesnt matter anyways)
        % if this cell is not a pyramidal cell, nan it out
        if ~contains(SuperRat(ses).units(i).type,'pyr')
            thispref=nan(size(thispref,1),size(thispref,2));
        end
        odorpref=[odorpref; thispref(bestind,:)];
    end
    %}
    % if we have two or more significant units
    if nansum(odorpref(:,3))>1
        % the cell indices for those who are selective
        ripinds=odorpref(odorpref(:,3)==1,4);
        odorinds=find(odorpref(:,3)==1);
        for j=1:length(ripinds)-1
            for k=j+1:length(ripinds)
                odorpairs(it,1)=ripcorrs(ripinds(j),ripinds(k));
                % absolute sum of the two, theyre signed so 0 is completely
                % opposite pref, and 2 is completely same pref, 1 is like
                % no relationship
                odorpairs(it,2)=abs(odorpref(odorinds(j,1))+odorpref(odorinds(k,1)));
                %odorpairs(it,2)=odorpref(odorinds(j,1))==odorpref(odorinds(k,1));
                it=it+1;
            end
        end
    end
    fprintf('Ses %d done \n',ses);
end

plot(odorpairs(:,1),odorpairs(:,2),'ko');
xlabel('Ripple Correlation');
ylabel('Odor selectivity correlation');
[~,m,b]=regression(odorpairs(:,1)',odorpairs(:,2)');
[h,p]=corr(odorpairs(:,1),odorpairs(:,2));
hold on;
plot([-.2 1],[-.2*m+b m+b],'r');
title(sprintf('%d pairs R^2 = %.2f ranksum p = %.4f',it,h,p));

%% okay so it looks like there is a breakout set of pairs who fire together

% given there are a minority of pairs who fire together, I guess one
% question is w


























%% heres the math behind the poisson decode

spkmat=nan(length(tstarts),length(unitdata),length(eventstarts));
h=waitbar(0,'gathering spikes');
for un=1:length(unitdata)
    allspikes=unitdata(un).ts;
    [~,trspikes]=event_spikes(allspikes,tstampmat(:),0,binsz);
    spkmat(:,un,:)=reshape(trspikes.*binsz,length(tstarts),length(eventstarts));
    % now fill each trialbin
    waitbar(un/length(unitdata),h,sprintf('gathering spikes for un %d',un));
end
close(h)
% now break up spkmat into three brain regions
% this is the expected firing rates on a given trial (mean across trials * binsize)
% 3d up this so we can do matrix math here
% 5hgis will be a nunits x ntsteps x ntsteps matrix

%spkpriors=nanmean(spkmat,3)'*.02; % multiply by timestep because the above fx divides
% now build the posterior distribution for each unit and each trial
for trial=1:length(eventstarts)
    % build from all other trials
    % so this is a 2d matrix of the mean firing rate for each cell at each
    % timepoint
    spkpriors=nanmean(spkmat(:,:,(1:length(eventstarts)~=trial)),3)'*.02; % multiply by timestep because the above fx divides
    
    % this trials rasters across units
    newevents=eventstarts(trial)+tstarts;
    % get spike counts for each bin
    trcounts=eventspikematrix(newevents',unitdata,0,binsz)'*binsz;
    % we want a square mat for each unit (p all times x each time
    % this trial)
    probmatL=nan(length(unitdata),length(tstarts),length(tstarts));
    decodedL=nan(length(tstarts),length(tstarts));
    % we'll loop through the templates
    
    h=waitbar(0,'starting');
    % for each timebin of this real trial data
    for timebin=1: length(tstarts) % for each timebin
        
        % probmat will be a m trial, x n unit x t timebin
        % p(t|s)=c* prod across units(  p(time)* fi(time)^ni / fact(ni)  * exp(Fi(time)))
        % calculate the probability of each time given this activity
        probmatL(:,:,timebin)=(spkpriors.^repmat(trcounts(:,timebin),1,size(trcounts,2)))./...
            repmat(factorial(trcounts(:,timebin)),1,size(trcounts,2)).*exp(-spkpriors);
        % now normalize so that probs for each unit sum to 1,
        % and then get the product?
        probmatL(:,:,timebin)=probmatL(:,:,timebin)./sum(probmatL(:,:,timebin),2); % make each units sum to 1
        % and line plot is what we ought to save for our priors
        decodedL(timebin,:)=prod(probmatL(:,:,timebin));% product across units
        
        % you can also log the theta angle on this trial too if you
        % want to plot the error of the decoder based on phase (it
        % hsould go from - to +
        waitbar(timebin/length(tstarts),h,'working');
    end
    close(h);
    
    decoded2=decodedL./nansum(decodedL);
    figure; imagesc(log(decoded2)'); set(gca,'CLim',[-10 -6]);
    
    % so we accumulate this across units, normalize somehow for
    % each unit and then get the decoded position for that trial
    
end

