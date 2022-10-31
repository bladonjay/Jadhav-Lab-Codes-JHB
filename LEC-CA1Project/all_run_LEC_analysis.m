%% LEC_CA1_DelayTask
%{
this is the general form the the lec ca1 paper



Figure 1 is on behavior exclusively, so we have distributions of events
locked to the delay onset.  This gives us a good idea of how long
everythign takes, and also whether there is any difference in decision
time or study time for correct vs incorrect trials

Figure 2 is histology, and a very quick PYR vs IN and CA1 vs PFC on burst
and acgram- gist is that ca1 is theta and LEC is not.  also firing rate
distributions overall

Figure 3 is object coding.  First some single cell examples of study and
test object coding, a quick plot of no object coding during delay, and then
cell counts, % study coders, % test coders, % time cells, % overlaps, maybe
a waterfall plot for study objects 4x4 showing normalized rates(I think
they will be normalized to the peak object, and be waterfalls sorted by the
peak, also one for the test objects too, could add those as waterfalls too
so its a 6x6, but the grids may be odd.

Figure 4 shoudl be a decoder, I think this will be interesting to dig into
incorrect trials and on timing, break down first and second test object
samples, and obviously show no coding during delay for either study or test
objects.  

figure 5 is all LFP, find where power is higher compared to pre study (take
last second or two of study, and compare to second or two before stuy
starts), then look at phase amp coh across and within regions.... if
anything pops out, look for spike phase coherence and overlap with coding
cells or INs

Figure 6, spike phase interactions if there are any, will have ot think
about that

Figure 7, swr figure if time, or its a cross-region spike field or coding
similarity analysis- see if coding peaks similarly across time or trials
for both regions the same, prob using the decoder

%}

% to load the file...
load('D:\Old Data\LEC Delay data\SuperRat5-29-20ALLLFP.mat'); % has ll th eLFP fiples
cd ('E:\GithubCodeRepositories\Jadhav-Lab-Codes\LEC-CA1Project');
%% a little about how the data were generated:
% 

edit run_behavior_stats