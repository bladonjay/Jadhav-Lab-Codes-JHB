%single trial coherence spectrum
function cs_singleTrialCoh(animal,day,epoch,trial,regions)
%regions should be 1x2 cell array with region strings
%CS33, day 3 is pretty good. trials 47, 31
%% Set Directories
[topDir, figDir] = cs_setPaths();
animDir = [topDir, animal,'Expt\',animal,'_direct\'];
%% Params
movingwin = [1000 30]/1000; 
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 40];
params.tapers = [3 2]; %DO NOT GO LOWER THAN THIS- MESSES UP 

%% Load Data
nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
tetinfo = loaddatastruct(animDir, animal, 'tetinfo');
rewardTimes = loaddatastruct(animDir, animal, 'rewardTimes', day);
coherence = loaddatastruct(animDir, animal,'coherence',day);

%Use all tetrodes
filt = ['isequal($area,''',regions{1},''') && $numcells > 0'];
tets1 = evaluatefilter(tetinfo{day}{epoch},filt);

filt = ['isequal($area,''',regions{2},''') && $numcells > 0'];
tets2 = evaluatefilter(tetinfo{day}{epoch},filt);

if isempty(tets1) || isempty(tets2)
    error('No tets found')
end

%% Get tet pair combinations
[A,B] = meshgrid(tets1,tets2);
pairs=reshape(cat(2,A',B'),[],2);
        
%% Loop over tet pairs
Coh = [];
for p = 1:size(pairs,1)
    tet1 = pairs(p,1);
    tet2 = pairs(p,2);
    
    disp(['Doing ', num2str(tet1), '-',num2str(tet2)]);
eegdat1 = loadeegstruct(animDir, animal, 'eeg',day,epoch,tet1);
eegdat2 = loadeegstruct(animDir, animal, 'eeg',day,epoch,tet2);

%% Get eeg data and time
betatime = geteegtimes(eegdat1{day}{epoch}{tet1});
eeg1 = eegdat1{day}{epoch}{tet1}.data;
eeg2 = eegdat2{day}{epoch}{tet2}.data;

%% Get trial time
trialstart = nosepokeWindow{day}{epoch}(trial,1)- 1;
npstart = nosepokeWindow{day}{epoch}(trial,1);
npend = nosepokeWindow{day}{epoch}(trial,2);
rewardstart = rewardTimes{day}{epoch}.allTriggers(trial);
trialend = rewardstart + 2; %do 1 second of reward? may need to adjust.

trialwin = [trialstart, trialend];

%% Use only eeg data within trial time
eegind = isExcluded(betatime,trialwin);
eeg1 = eeg1(eegind);
eeg2 = eeg2(eegind);

%% Calc coherence
[coh,Phi,~,~,~,t,freq] = cohgramc(eeg1,eeg2,movingwin,params); 
coh = coh';

Coh = cat(3,Coh, coh);
end

Coh = mean(Coh,3);

%% Zscore
mn = coherence{day}{epoch}.mean;
sd = coherence{day}{epoch}.sd;

%Coh = bsxfun(@rdivide,(Coh - mn),sd);
%% Plot
close all
std = 5;
    s = gaussian2(std,(3*std));
    cohgram = filter2(s,coh,'valid');
    figure, 
    colormap(jet);
    imagesc([(trialstart-npstart):(trialend-npstart)], [1:max(freq)],cohgram)
    set(gca,'YDir','normal')
    colorbar
    hold on
    plot([0 0],[1 max(freq)],'k--', 'LineWidth', 1.5);
    plot([npend-npstart npend-npstart], [1 max(freq)],'k--', 'LineWidth', 1.5);
    plot([rewardstart-npstart rewardstart-npstart], [1 max(freq)],'k--', 'LineWidth', 1.5);
    
    
