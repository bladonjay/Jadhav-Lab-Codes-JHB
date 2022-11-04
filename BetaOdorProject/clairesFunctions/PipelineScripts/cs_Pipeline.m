%extract data
%matclust
%position track
%preprocess
    %create pos, lfp, spikes, dio
    %create task, tetinfo, cellinfo and update all with tags
    %extract ripples
    %lfp filters
    %odorTriggers files
%spectrograms
    %baseline specgrams
    %specgrams
%coherence, cohereograms
%cell type tag
%odor selectivity tag
%odor psths
%phase locking tag
clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

edit preprocess_CS31
preprocess_CS33
preprocess_CS34
preprocess_CS35
preprocess_CS39
preprocess_CS41
preprocess_CS42
preprocess_CS44

%
%% --- Behavior --- %%
cs_getNosepokeWindow(animals)

cs_performanceScatter(animals);

cs_turnFromNPDirection
%% --- Spectrograms --- %%
DFScs_eventTrigSpecgram
% notes on this function:
% core function is mtspecgramc with inputs:
% first is full day lfp
% params: tapers=[1 1]; movingwin = [1000 20]/1000 fs=1500

%% --- Coherograms --- %%
DFScs_eventcoherence

%% --- Cell Data --- %%
edit cs_listAllCells; % i think it gets you all cells

edit cs_cellTypeTag; % pyr or in

edit cs_listRunCells;

% old np and odor responsive fx, defunct
% cs_listNPcells
% cs_odorResponsiveCells



% this is the current version we are using.  This will list np cells and
% those that fire at least once per trial during odor period
edit cs_listNPCells_v2; % gets task Responsive cells (this now includes odor responsive)
%%
%{
ran npcells first
cs_cellselectivitytag: 49 & 53
run
cs_listselectivity: 49 & 53

when you run selectivitytag you get more cells
but for listselectivecells, that list has overlap with npcells has 44 and
52


%}

%% --- Odor Selectivity --- %% 
% this gathers all the selective cells
edit cs_cellSelectivityTag;
edit cs_intSelectivityTag;
% this plots all of them

% these plot the red to blue plots and the xcorrs correct incorrect
edit cs_odorSelectivity_v2
edit cs_odorSelectivity_INs


%% plot some cells
% this plots all cells
edit cs_plotRasterPSTH_v3

% this plots only odor selective cells
edit cs_plotRasterPSTH_specifiedCells_v2


%% --- PCA/ Population Analysis ---%
cs_PDI(win,binsize);

cs_individualTrialPDI
cs_PDIBehaviorCorrelation

cs_PCA(region, win, binsize, selectiveonly);
cs_PCA_3d(region, win, binsize, selectiveonly);

%% --- Phase locking --- %%

% here is where Jay has updated most of these functions

%edit jhb_calcBetaRRallTets % this gets the beta and rr filtered LFPs, unnecessary to run though

edit cs_phaseLocking % this version has no 'high beta' windows
edit cs_phaseLocking_IN

edit cs_listPhaseLockedCells %
edit cs_listPhaseLockedINs.m

%%
% then all cells
edit cs_cellPopulations
%%

%cs_phaseLocking_v3

edit cs_getHighBetaTimes % also unnecessary, we will not be using high beta as a filter
edit cs_phaseLocking_highBeta

% this generates population figures, gets dial plots for all pfc, all ca1,
% then proportion of units that pass rtest for each region, then the dial
% plots for within and cross coherence
edit cs_popPhaseLock_v2
edit cs_popPhaseLockCells_v2

edit cs_phaseLocking_performance_v3
edit cs_phaseLocking_performance_INT

%%
% plot the phase locking and the interaction between selectivity and
% locking
jhb_plotAllProportions

