%% 1.
% Figure out the behavioral parameters during odor presentation and 
% recall. From average movement/speed curves (similar to code I had sent), 
% define a time point beyond which the animal had definitely moved out of 
% the nose-poke. We cannot assign any neural correlates to odor or choice 
% selectivity beyond this point. 

%CS- nosepoke duration based on DIOs- find when nosepoke went into off state,
%first instance after established odor trigger times. Then can calc avg for
%animal, and across animals. 
animals = {'CS31','CS33','CS34','CS35','CS39'};
cs_getOdorTriggers_v2(animals)
meanNPDuration = cs_meanNosepokeDuration(animals);
mean(meanNPDuration)

%% 2.
%Once you define the "average" time-point 
% based on all trials, also notice the variance of the motion curves. A 
% more conservative criterion can be [average - 1 standard deviation].  
% Alternatively, for each trial, it will be useful it will be useful to 
% define a time-point when animal moved away from nose-poke using the 
% motion curves from just that trial, since this will be variable across 
% trials. I expect we will see selective/preparatory activity that will 
% peak just before the animal moves out, and it will be powerful if we can 
% show correlation between these two measures.
% 
% Do a similar analysis for population activity of selective cells, for 
% CA1 and PFC spearately. You are defining a time point of separation of 
% odor1 vs odor 2 curves with average across all trials. Now for each 
% individual correct trial, you know that the stimulus is Odor 1 or Odor 
% 2. Can you find out the time point when population activity in the 
% current trial hits the average (or average+1 std) selectivity point of 
% Odor 1? Meaning: you know that on average, the significant separation of 
% curves happens at x, say 400ms. For the metric that you are using for 
% this average analysis, either Euclidean distance or PCs, you know the 
% value of the average curve for Odor 1 and Odor 2 at 400ms - they will 
% separated by some distance since they are significantly different. Can 
% you find out when the population activity for individual trials hits 
% that average value? It is arguable that this is the time for neural 
% separation in individual trials, and we can see how it correlates to the 
% animal's behavioral moving out time. I expect the neural separation will 
% be >100ms than behavioral separation on a trail-by-trial basis. We may 
% also see some CA1 vs PFC differences here that are not apparent in 
% averages. If we see a reasinable result, we can repeat this for error 
% trials as well.
%
%CS - compare the PV from one trial to the average across all trials for
%the opposite odor (from the same animal...?). Calc the PDI between these
%two vectors. At what time point does the distance between them hit the
%minimum significant distance between the average vectors across
%trials/animals? Do for all individual trials, then compare each of these
%timepoints to the timepoint at which the animal left the nosepoke.

cs_individualTrialPDI([0.2 1.5],0.1);
cs_individualTrialPDI_v2([0.2 1.5],0.1);
cs_PDIBehaviorCorrelation;

%% 3.
% Do a similar analysis for time-points for single cell selectivity. For 
% each selective cells, can you define a time point for when separation of 
% curves is significant (use smallest time bins that you can use)? Get an 
% average/distribution of this time point across cells, separately for CA1 
% vs PFC.

cs_singleCellSelectivityOverTime(win, binsize)
cs_cellSelectivityEmergenceDistribution

%%
% Do a cross-correlation (actually cross-covariance) for CA1 vs PFC 
% neurons during the recall period. We have done this, corss-covariance 
% during theta, in my Neuron 20016 papers and also the J. Neursoci paper 
% last year. You can do a similar analysis during just the recall period, 
% and we can see if there is any directionality in the peak of the 
% corss-covariance. Eventually, we can supplement the directionality 
% analysis using GLM, Granger causality, or a Information theoretic metric.
cs_xcorr

%try creating a histogram of where the peak is for each pair, rather than
%averaging over all
cs_xcorr_v2

% We can leave the jPCA analysis for later.
% 
% 
% Switching to LFP:
% 
% - Look at trial-by-trial LFP (for each trial, you can plot in a single 
% figure, raw traces, filtered in beta band, filtered in theta/HRR band), 
% and visually inspect for any trends - a) do beta oscillations look 
% variable in onset after nosepoke, b) any clue for respiration frequency 
% in raw LFP? We can expect that beta oscillations will start after 
% sniffing onset, so check if there is simultaneous increase in beta and 
% some other frequency, possibly theta
cs_trialLFP;
% 
% - Plot average spectrogram for lower band (~0-20Hz) with better 
% frequency resolution, and check if any sniffing band is apparent. OB 
% would be most promising here.
% 
% - Are trial-by-trial spectrograms too noisy, or can we find any sniffing 
% band in these plots that may be variable across plots. OB would be most 
% promising here.
% 
% - Phase-locking. If beta onset is very variable across trials, 
% phase-locking data should be limited to periods after beta onset till 
% movement out of nose-poke (high beta power filter)
cs_phaseLocking_highBeta;

% 
% - What beta phases do CA1 and PFC neurons prefer, using this more 
% restricted high beta power filter? We can also do a cross-covariaince 
% analysis as above (for spikes) using this high beta power filter. I 
% suspect we can either use a high beta power filter, or take all time 
% periods after a "beta onset" time, if we can define such a time for each 
% trial.
cs_popPhaseLock_v2
% 
% - Checking for phase-locking properties of single cells for a possible 
% sniff/theta frequency may not be a bad idea.
% 
% - The jPCA analysis will allow us to separate the oscillatory population 
% activity from non-osc. So it may be possible to check selectivity of 
% just the osc part. But lets find the LFP correlate of the sniff 
% frequency first! If I remember correctly, Igarashi et al, did not deal 
% with this either. It may be possible that in a discrimination task like 
% this, it is not a prominent component of CA1 or PFC LFP. But it has to 
% be there in the OB LFP.


% Selectivity Emergence vs NP exit on error trials
% Trial by trial selectivity index
% Session by session PV divergence/NP offset
% Refine PV divergence/NP offset 
    %spiking criterion for each trial, 
    %look at good vs bad correlation trials and figure out what the
    %difference is, can further refine from there
    %also compare to error trials
    %also add beta onset, or refine to beta onset
%Place fields of selective neurons 
    %plot trajectories independently? 
%Wavelet spectrograms
%Cell phase preference plots (rose plots?)
%Xcorr for only phase locked neurons, limit time window to <100ms? 
%area under curve to determine if there is a bias