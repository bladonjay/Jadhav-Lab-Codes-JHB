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

preprocess_CS31
preprocess_CS33
preprocess_CS34
preprocess_CS35
preprocess_CS39
preprocess_CS41
preprocess_CS42
preprocess_CS44
%% --- Behavior --- %%
cs_getNosepokeWindow(animals)

cs_performanceScatter(animals);

cs_turnFromNPDirection
%% --- Spectrograms --- %%
DFScs_eventTrigSpecgram

%% --- Coherograms --- %%
DFScs_eventcoherence

%% --- Cell Data --- %%
cs_listAllCells;

cs_cellTypeTag;

cs_listNPCells;


%% --- Odor Selectivity --- %% 
cs_cellSelectivityTag;

cs_odorSelectivity_v2;

cs_plotRasterPSTH_v3
cs_plotRasterPSTH_specifiedCells_v2


%% --- PCA/ Population Analysis ---%
cs_PDI(win,binsize);

cs_individualTrialPDI
cs_PDIBehaviorCorrelation

cs_PCA(region, win, binsize, selectiveonly);
cs_PCA_3d(region, win, binsize, selectiveonly);

%% --- Phase locking --- %%
%cs_phaseLocking_v3
cs_getHighBetaTimes
cs_phaseLocking_highBeta

cs_popPhaseLock_v2
cs_popPhaseLockCells_v2

cs_selectivityPhaseLockOverlap_v3

cs_phaseLocking_performance_v2