%% ClaireWrapper

% first build your dataset;
%load('E:\Claire Data\ClaireData22-Apr-2020.mat')
%% or to rebuild the dataset from scratch
edit ClaireDataAggregator

%% or load the dataset like this:
% this dataset contains 


% now you can linearize the positional data
% as currently is 2-10-2020 this will ask you to redraw cs39 over again
% this is because this is a partial long track, remember no 90* angles on
% your tracks!
edit LinearizePositionTmaze

% the 4 traj is better if your runs are all wonky (like for short tracks)
edit LinearizePositionTmazeNTraj
% some helpful fx for that:
 
% - one is getcoord_Tmaze, this lets you draw the trajectories, or at
% least is a skeleton for opening a gui and drawing some mean trajectories
%%  if you've already made linearized position data

%load('E:\Claire Data\ClaireDataFull-LinPos-LongTrack.mat');
% removes short sessions
SuperRat(~logical([SuperRat.longTrack]))=[];


%%
load('E:\Claire Data\ClaireDataFull-LinPos-LongTrack-ObjCoding.mat')

% first parse pyrams and interneurons based on rate and burst index
edit ParsePyrIN

% to plot some place fields out use this:
edit ClaireQuickPlacePlot
% plot linearized place fields
edit PlotLinTrajectories
% calculate them for the ripple decoding and to crossref with obj coding
edit CalcLinTrajectories %
% get object coding
edit Object_selectivity

    
% get spatial info for object cells
edit PlotObjPlaceCells

%% and the meat of the figures
edit SummaryObjPlaceInteractions

edit Claire_Odor_Routedecoder
%%  This is the beta coherence datast

load('G:\Claire Data\ClaireData22-Apr-2020.mat')
betaDir='E:\ClaireData';

edit ClaireBetaAnalysis

% 
% rr analysis just has 
edit ClaireRRAnalysis

edit ClaireCrossCoherence

%% now to analyze these data
% one fun question is whether odor representations are replayed with goal
% location representations, although it may be tough as there are only four
% options, so the basrates are very high
edit ClaireRippleAnalysis

%%
clearvars -except SuperRat



%% More notes
%{
notes for splitter:
wenbo did a correlation type scoring across the two linearize trajectories
One thing yu can do is run an xcorr between the two rate maps.  The other
way to do this is to bootstrap the absoluite differenes between the rate
maps. the algorithm would go as follows:"
- gather the full ratemaps for the real data, and at each bin take the
absolute difference between the two maps
- to get signficance, gather each 'run' and randomize which trajectory it is.
Then you can build rate maps for each trajectory, and then run the absolute
difference between the two.  WEnbo ran the diff over sum, but
mathematically, that underrepresents high firing pixels and i dont know if
i want that...
the regression i think we're looking for is how similar the rate maps are,
but i'll need to normalize them by maybe their rates overall?


%}

%{
Notes from 4/10/20
for figure 5- put letters for subpanels
-remove odor 1 vs odor 2 rates for cells that are only responsive and not
selective
-add raster for odor responsiveness for some cells
-more example cells
-in supplement, show specifically cells that are odor selective and show
their route selectivity is +, mnopt selective, or opposite selecting
-for pfc, cause the regression works, pick cells that have matching route
selectivity and odor selectivity
-for route selectivity, nspikes per run is fine
-for pfc place fields need to be called spatial fields
-use a binomial test to compare the probability of responsive
cells/selective cells having pfs to those not responsive or selective
-replace the bar graph with pf numbers as x axis with the regression plot
of rates during run and rates during odor presentation
-use claires glm code to decode routes given odor activity.


%}


%{
glm code
% ok so Cellmatrix is n events by m cells, and its the number of spikes per
event row.  my window will be 800 msec.  So the goal here will be to use
exactly her decoder and then apply it to multiple spatial bin sizes.

K = 5;
cv = cvpartition(numtrials, 'kfold',K);
%        
for k=1:K
    % training/testing indices for this fold
    trainIdx = cv.training(k);
    testIdx = cv.test(k);
    
    % train GLM model
    
    warning('off','all');
    %disp('Fitting GLM')
    mdl = fitglm(Cellmatrix(trainIdx,:), trialIDs(trainIdx),'Distribution', 'binomial');
    
    % predict regression output
    Y_hat = predict(mdl, Cellmatrix(testIdx,:));
    
pred = round(Y_hat); %glm outputs very small values instead of zeros (binomial) for some reason, round to zero
ids = trialIDs(testIdx);
fract_correct(k) = sum(pred == ids)/length(pred);

end
fract_correct = mean(fract_correct);

%Do shuffling
for s = 1:100
    if mod(s,10) == 0 || s == 1
    %disp(['Doing shuffle- iteration ', num2str(s)])
    end
    
    shuff_trialIDs = trialIDs(randperm(length(trialIDs)));
    mdl = fitglm(Cellmatrix(trainIdx,:), shuff_trialIDs(trainIdx),'Distribution', 'binomial');
    Y_hat_shuff = predict(mdl, Cellmatrix(testIdx,:));
    pred_shuff = round(Y_hat_shuff);
    ids_shuff = shuff_trialIDs(testIdx);
    
    fract_correct_shuff(s) = sum(pred_shuff == ids_shuff)/length(pred_shuff);
end
fract_correct_shuff = mean(fract_correct_shuff);

%}
%% random sketchpad

odorinfo=[SuperRat(ses).trialdata.sniffstart SuperRat(ses).trialdata.leftright10,...
    SuperRat(ses).trialdata.CorrIncorr10];
rundata=SuperRat(ses).LinCoords;
% now lets go trial by trial and see what the odor ID is, its correct
% incorrect, and the route the animal took.

% Correct: l is L, 0 is R
%
figure; 
for i=1:5
subplot(5,1,i);
trial=i;
runstart=find(rundata(:,1)>odorinfo(trial,1),1,'first');
runstop=runstart+find(rundata(runstart:end,8)>55,1,'first');
plot(rundata(runstart:runstop,2),rundata(runstart:runstop,3));
title(sprintf('odorID=%d, corrincorr=%d',odorinfo(trial,2),odorinfo(trial,3)));
end

