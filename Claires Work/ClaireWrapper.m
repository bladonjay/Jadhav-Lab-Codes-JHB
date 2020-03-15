%% ClaireWrapper

% first build your dataset;
load('E:\Claire Data\ClaireData 3-4-2020.mat')

edit ClaireDataAggregator

%% or load the dataset like this:
% this dataset contains 
load('E:\Claire Data\ClaireData 2-12-2020.mat')

% now you can linearize the positional data
% as currently is 2-10-2020 this will ask you tyo redraw cs39 over again
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
edit CalcLinTrajectories
% get object coding
edit Object_selectivity

    
% get spatial info for object cells
edit PlotObjPlaceCells

edit SummaryObjPlaceInteractions

%% now to analyze these data

edit ClaireRippleAnalysis

%%
clearvars -except SuperRat
