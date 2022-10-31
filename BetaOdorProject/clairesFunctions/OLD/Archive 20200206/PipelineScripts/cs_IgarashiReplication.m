%IGARASHI ET AL REPLICATION -- MASTER SCRIPT

%% --- Global Params ---%%
[topDir, figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC','OB'};

%% --- Figure 1d, Average Speed ---%%
% Get average trial length
avgTrialLength = cs_avgTrialLength(animals, topDir);
%avgTrialLength = round(avgTrialLength); 

% Get average speed and plot
cs_avgTrialSpeed(animals, topDir, avgTrialLength, figDir);

%% --- Figure 1e, Normalized Power Spectra --- %%
% Percent change from baseline (0-150ms before nosepoke)
DFScs_eventTrigSpecgram_percentChange
cs_plotSpecgram('eventTrigSpecgramData_CA1_percentChange_low_-1000-2000ms_tapers1-1_20msBins.mat')
cs_plotSpecgram('eventTrigSpecgramData_OB_percentChange_low_-1000-2000ms_tapers1-1_20msBins.mat')
cs_plotSpecgram('eventTrigSpecgramData_PFC_percentChange_low_-1000-2000ms_tapers1-1_20msBins.mat')

%% --- Figure 1f, Positions with High Beta Power --- %%

%get average beta power across tetrodes for each day and region
cs_avgBetaPower(animals, topDir, regions)

%Create spatial map with beta power
cs_betaPowerPos(animals, topDir, regions, figDir)

cs_betaPowerPos_nponly(animals, topDir, regions, figDir, 1);
cs_betaPowerPos_nponly(animals, topDir, regions, figDir, 0);


%detect beta oscillatory episodes
regions = {'CA1','PFC','OB'};
cs_detectOscillatoryEpisodes(topDir, animals, regions)

%Beta Coherene
DFScs_epochCoherence

regions = {'PFC-OB','CA1-OB','CA1-PFC'};
cs_coherencePos(animals, topDir, regions, figDir);
cs_coherencePos_nponly(animals, topDir, regions, figDir, 1);
cs_coherencePos_nponly(animals, topDir, regions, figDir, 0);
cs_coherencePos_nponly(animals, topDir, regions, figDir, 2);
%% --- Figure 1i, Coherence spectrum at 200ms after nosepoke as function of Frequency ---%%
cs_instantaneousCoherence(animals, topDir, regions, figDir);

%% --- Figure 1j, percentage of high coherence times as function of trial time ---%%
cs_percentHighCohTrialTime(topDir, animals, regions, figDir)

%% --- Figure 2b, coherence spectrum as function of frequency for correct vs incorrect trials ---%% 
cs_instantaneousCoherence_CI(animals, topDir, regions, figDir)

%% --- Figure 2d, cross-frequency coherence --- %%



%% --- Figure 4a, Odor Specific Representations --- %%
regions = {'CA1','PFC'};

cs_cellTypeTag(topDir, animals, regions)

win = [0.6 0.6];
binsize = 0.1;

cs_odorSelectivity_v2(topDir, figDir, animals, regions, win,binsize)

%% --- Figure 4c,  PDI --- %%
win = [0 1.5];
binsize = 0.05;

cs_PDI(win,binsize)


