% OUTLINE
%{
extract data
matclust
position track
preprocess
    create pos, lfp, spikes, dio
    create task, tetinfo, cellinfo and update all with tags
    extract ripples
    lfp filters
    odorTriggers files
spectrograms
    baseline specgrams
    specgrams
coherence, cohereograms
cell type tag
odor selectivity tag
odor psths
phase locking tag
%}
%% must be in the directory that this script is in
clear

% first add the repositories above this directory EXCEPT FOR 'OLD'
% will need to add 'general repo' as well...
jhb_setBetaOderpaths

% now set our data path
cs_setPaths



%% this is generated pretty early and has a lot of hardcoded data
% given the data repository in its current form this is unnecessary

%{
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
%}

% the above requires the following files:
% A directory with the 'topdir' then animalExpt, then animal_direct
% example:
% E:\BetaOdorProject\CS31Expt\CS31_direct
% with files:
% 

% the above generates the following files:

%% --- Behavior --- %%
% all 8 animals
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

edit cs_getNosepokeWindow
%cs_getNosepokeWindow(animals)

edit cs_performanceScatter;
%% and get a turn direction

edit cs_turnFromNPDirection
%% --- Spectrograms --- %%

% THIS FUNCTION DOES NOT CURRENTLY RUN, I AM WORKIGN ON ADDING SOURCECODE
% AND LFP DATASET

% this function is not working, some source code is outdated specifically
% wrt filterframework files, I ended at animaldef in the src folder, which
% has shantanus old code in it, not claires new animal defs'
edit DFScs_eventTrigSpecgram
% notes on this function:
% core function is mtspecgramc with inputs:
% first is full day lfp
% params: tapers=[1 1]; movingwin = [1000 20]/1000 fs=1500



%% --- Coherograms --- %%
% THIS ALSO DOESNT RUN, DO NOT TRY....
edit DFScs_eventcoherence

%% --- Cell Data --- %%
edit cs_listAllCells; % i think it gets you all cells

edit cs_getSpikeWidth;

edit cs_cellTypeTag; % pyr or in

edit cs_listRunCells; % lists cells that fire during run epochs

% old np and odor responsive fx, defunct
% cs_listNPcells
% cs_odorResponsiveCells

% this is the current version we are using.  This will list np cells and
% those that fire at least once per trial during odor period
edit cs_listNPCells_v2; % gets task Responsive cells (this now includes odor firing)
edit cs_listNPInt; % gets task responsive and now odor firing
%% notes on discrepancies between versions:
%{
ran npcells first, she updated np cells after runing the rest of the code
though, so we need to rerun here down after the correct npcells is run
cs_cellselectivitytag: 49 & 53
run
cs_listselectivity: 49 & 53

when you run selectivitytag you get more cells
but for listselectivecells, that list has overlap with npcells has 44 and
52


%}

%% --- Odor Selectivity --- %% 
% this calculates all the selective cells
edit cs_cellSelectivityTag; % tags selective cells, contains listSelectiveCells
edit cs_intSelectivityTag; % tags selective ints, contains listSelectiveInt
% this plots all of them

% these plot the red to blue plots and the xcorrs correct incorrect
% need to put the regression curves there (linear model)
edit cs_odorSelectivity_v2
edit cs_odorSelectivity_INs


%% plot some cells
% this plots all cells
edit cs_plotRasterPSTH_v3

% this plots only odor selective cells
edit cs_plotRasterPSTH_specifiedCells_v2


%% --- PCA/ Population Analysis ---%

% for current paper, run this function and not what's below it:

edit Run_PCA_SJ_copy.m


% havent vetted this yet
cs_PDI

cs_individualTrialPDI
cs_PDIBehaviorCorrelation

win = [.2 1];
binsize = 0.1;
selectiveonly=0;
regions={'CA1','PFC'};
for i=1:2
    region=regions{i};

    cs_PCA(region, win, binsize, selectiveonly);
    cs_PCA_3d(region, win, binsize, selectiveonly);
end
%% --- Phase locking --- %%


% this gets the beta and rr filtered LFPs, unnecessary to run though
% edit jhb_calcBetaRRallTets 

% calculates phase locking for each cell
edit cs_phaseLocking % uses current np cells
edit cs_phaseLocking_IN % uses np cells old

% used for later, builds reference lists
edit cs_listPhaseLockedCells % generates a list in acrossanimals analyiss
edit cs_listPhaseLockedINs.m% same as above


%%
% then all cells, this generates a struct containing lists of cells for selectivity,
% responsiveness and coherence 
edit cs_cellPopulations
%%

% she went down the 'find high beta' rabbit hole, dont use
%edit cs_getHighBetaTimes % also unnecessary, we will not be using high beta as a filter
%edit cs_phaseLocking_highBeta



% this generates population figures, gets dial plots for all pfc, all ca1,
% then proportion of units that pass rtest for each region, then the dial
% plots for within and cross coherence
edit cs_popPhaseLock_v2 % this version works with dial plots
edit cs_popPhaseLockInt_v2 % current version, dial plots (and current task responsive


% these two work now as well (there are other versions for all cells not
% just coherent ones...
edit cs_phaseLocking_performance_v4
edit cs_phaseLocking_performance_INT

%%
% plot the phase locking and the interaction between selectivity and
% locking
edit jhb_plotAllProportions

