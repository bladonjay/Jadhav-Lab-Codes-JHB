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

preprocess_CS31
preprocess_CS33
preprocess_CS34
preprocess_CS35


%% --- Behavior --- %%
cs_getOdorTriggers_ss(prefix, days,epochs)

cs_stateSpacePerf_plotAllAnimals(animals, daycell) %bar graph of performance
cs_stateSpacePerf_plotAllAnimals({'CS31','CS33','CS34','CS35','CS39'}, {[1,2,3],[1,2,3,4],[1,2,3,4],[1,2,3]})
%% --- Spectrograms --- %%

DFScs_eventTrigSpecgram
cs_plotSpecgram(eventTrigSpecgram_file, varargin)

%% --- Coherograms --- %%
DFScs_eventcoherence
cs_plotCoherence(cohDataFile)

%% --- Cell Data --- %%
cs_listAllCells;

cs_cellTypeTag;

cs_listNPCells;


%% --- Odor Selectivity --- %% 
cs_cellSelectivityTag;
cs_listSelectiveCells;

cs_odorSelectivity_v2; % this makes the rasters and generates the anticorrelation

cs_plotRasterPSTH_v3
cs_plotRasterPSTH_specifiedCells_v2


%% --- PCA/ Population Analysis ---%
cs_PDI(win,binsize);

cs_PCA(region, win, binsize, selectiveonly);
cs_PCA_3d(region, win, binsize, selectiveonly);

%% --- Phase locking --- %%
%cs_phaseLocking_v3
cs_getHighBetaTimes
cs_phaseLocking_highBeta
%calls cswt_beta_phaselocking and cs_listPhaseLockedCells
cs_cellPhaseLockTag

cs_popPhaseLock_v2
cs_popPhaseLockCells_v2

cs_selectivityPhaseLockOverlap_v2