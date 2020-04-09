%% ClaireWrapper

% first build your dataset;
load('E:\Claire Data\ClaireData4-6-20.mat')

edit ClaireDataAggregator

%% or load the dataset like this:
% this dataset contains 
l

% now you can linearize the positional data
% as currently is 2-10-2020 this will ask you to redraw cs39 over again
% this is because this is a partial long track, remember no 90* angles on
% your tracks!
edit LinearizePositionTmaze

edit LinearizePositionTmaze4Traj
% some helpful fx for that:
 
% - one is getcoord_Tmaze, this lets you draw the trajectories, or at
% least is a skeleton for opening a gui and drawing some mean trajectories
%%  if you've already made linearized position data

%load('E:\Claire Data\ClaireDataFull-LinPos-LongTrack.mat');
% removes short sessions
SuperRat(~logical([SuperRat.longTrack]))=[];


%%
%load('E:\Claire Data\ClaireDataFull-LinPos-LongTrack-ObjCoding.mat')

% first parse pyrams and interneurons based on rate and burst index
edit ParsePyrIN

% to plot some place fields out use this:
edit ClaireQuickPlacePlot
% plot linearized place fields
edit PlotLinTrajectories
% calculate them for the ripple decoding and to crossref with obj coding
edit CalcLinTrajectories % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% get object coding
edit Object_selectivity

    
% get spatial info for object cells
edit PlotObjPlaceCells

edit SummaryObjPlaceInteractions

%% now to analyze these data

edit ClaireRippleAnalysis

%%
clearvars -except SuperRat


%%
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

