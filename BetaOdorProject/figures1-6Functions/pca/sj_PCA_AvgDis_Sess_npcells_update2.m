% SJ
% June 2021. Precursor to running from trialSpikes. Directly use npCells
% with Pyr and Interneurons, and get PC Discrimin per Sess/Day

% Can Remove Extraneous Stuff

% June 2021. Use the trialspikes structure from Claire, with both Pyr cells
% and Interneurons

% Similar to sj_PCA_AvgDis_Sess1

%edited 12/12/2019
%calculates PCA for left and right trials. Can either consider all trials
%individually, or just take mean firing rate of each cell across trials. 
%plots new vectors for PCs 1-3. Then performs a shuffle to determine at
%which bins the right and left vectors significantly diverge. 

%function sj_PCA1(region, win, binsize, selectiveonly)
% SJ: Make it a script

%close all
%sj_PCA('CA1',[-0.2 1], 0.1, 0)

%% --- Params 
% animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
% [dataDir, figDir] = cs_setPaths();
% figdir = [figDir,'PCA\']; 

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};


% for PFC example
regionIDX=2; animno = 8; day = 2;
% for HPC example
regionIDX=1; animno = 3; day = 3;


region = regions{regionIDX}; 
animal = animals{animno};
%dataDir = '/Users/Shantanu/data25/OLF_CS/Data/';
%figdir = '/Users/Shantanu/data25/OLF_CS/Data/sj_Figures';
dataDir='E:\Brandeis datasets\Claire Data\';
animalDir = [dataDir,animal,'Expt\'];
load([animalDir,animal,'cellinfo.mat'])
%daystr = getTwoDigitNumber(day);
daystr=sprintf('%02d',day);
load([animalDir,animal,'nosepokeWindow',daystr,'.mat'])

% jays loading
%load('E:\Brandeis datasets\Claire Data\CS33Expt\CS33nosepokeWindow01.mat')
%load('E:\Brandeis datasets\Claire Data\CS33Expt\CS33cellinfo.mat')
%

win = [0.2 0.9]; % Implies from -0.2 to 1
binsize=0.1;
binsize_ms=1000*binsize;
selectiveonly=0;
%timeaxis = -1000*win(1):binsize_ms:1000*win(2);

digits(16);

win1 = 0 - win(1); %make negative

winstr = [(num2str(win1*1000)),'-',num2str(win(2)*1000),'ms'];

if selectiveonly == 1
    selstr = '_selectivecells';
else
    selstr = '';
end

%% --- Calculate
disp('creating column vectors')
[leftfr, rightfr, lefttrials, righttrials, cellinds, np_left,np_right,np_all] = sj_columnVectors_day_npCells(animal, animno, day, region, win, binsize);
%[leftfr, rightfr, lefttrials, righttrials, cellinds] = cs_columnVectors(animals, region, win, binsize, selectiveonly);
%
% Get Bins
allbins = -0.5:binsize:1.5-binsize;
goodbins = win1:binsize:win(2)-binsize; 
%binind = lookup(goodbins,allbins);
%numtimebins = length(goodbins);

numtimebins = length(goodbins); 

%% 


% PCA on Firing rate 
% ----------------
LRcombined = [leftfr; rightfr];
[dim, num_data] = size(LRcombined);  
% dim/trials/obs/datapoints = Ntrials (left+right)
% num_data = number of neurons
% LRcombined or X is a Ntrials x Nneurons matrix

leftidx = 1:size(leftfr,1);
rightidx = size(leftfr,1)+1:size(LRcombined,1);

[COEFF, SCORES, LATENT, TSQUARED, EXPLAINED, MU] = pca(LRcombined);
% COEFF returns N_neurons prncipal components in columns, with each PC
% having Nneurons coefficients along original dimensions

% SCORES are the projections: Original data in Ntrials projected along
% the new N_neu dimensions - we will look at the first 3


%% 

% PCA on firing rate vectors
% ------------------------

leftvector=[]; rightvector=[]; %For each timebin, save NtrialXNnenuron firing matrix

Nneu = size(lefttrials,1);
Nlefttr = length(leftidx);
Nrighttr = length(rightidx);
Ntotaltr = Nlefttr + Nrighttr;

bins = -win(1):binsize:win(2);
Nbins = length(bins)-1; % histogram will have bins-1 edges

% For each bin, get a NtrxNneu matrix, and then do PCA

for i=1:Nbins
    
    storeleftbindata = []; storerightbindata = []; %Reinitialize store for current bin
    for neu = 1:Nneu
        currleft = lefttrials{neu};
        currright = righttrials{neu};
        
        currleftbindata = currleft(:,i); %current neurons data for timebin i
        storeleftbindata = [storeleftbindata, currleftbindata];
        
        currrightbindata = currright(:,i); %current neurons data for timebin i
        storerightbindata = [storerightbindata, currrightbindata];    
    end
    
    leftvector{i} = storeleftbindata;  % each bin has NtrXNneu data
    rightvector{i} = storerightbindata;
    
end


%% 



% a) Do PCA in each timebin
% ---------------------
% b) Also, GET average PCA trajectory for left and right trials, and also per trial
% Also get DISTANCE METRIC FOR 1st 3 PCs only

% For storing PC values
% -------------------
store_Scores = []; store_Explained = [];
store_Scores3 = []; store_Explained3 = [];

% For storing Average and variance of PC values
% ---------------------------------------------
meanleft=[]; meanright=[];
lefttrials_pc = zeros(Nlefttr,Nbins,3); righttrials_pc =zeros(Nrighttr,Nbins,3);

% Do only DIST between average trajectories for Left and Right
store_avg_DISTPC1=[];  % Avg left vs Avg right
store_avg_DISTPC=[];  % Avg left vs Avg right




for i = 1:Nbins
    X_left = leftvector{i};
    X_right = rightvector{i};
    
    LRcombined = [X_left; X_right];
    [COEFF, SCORES, LATENT, TSQUARED, EXPLAINED, MU] = pca(LRcombined);
    
    store_Scores{i} = SCORES;  % Scores = NtrXNneu
    store_Explained{i} = EXPLAINED;
    store_Explained3 = [store_Explained3;sum(EXPLAINED(1:3))]; % explain variance in each time bin
    store_Scores3{i} = SCORES(:,1:3); % First 3 components only. Ntrx3
    
     % Store for current bin - mean and err for PCs
     % --------------------------------------------
    curr_scores = store_Scores3{i};
    currleft_scores = curr_scores(leftidx,:);
    currright_scores = curr_scores(rightidx,:);
    meanleft(i,:) = nanmean(currleft_scores,1);   % For current timebin, save mean across trials of 1st 3 PCs/1st PC
    meanright(i,:) = nanmean(currright_scores,1); 
    errleft(i,:) = nansem(currleft_scores,1);     % For current timebin, save sem across trials of 1st 3 PCs/ 1st PC
    errright(i,:) = nansem(currright_scores,1);
    stdleft(i,:) = nanstd(currleft_scores,1);     % Same for st dev
    stdright(i,:) = nanstd(currright_scores,1);    
    
    
    % Gather to plot trial-by-trial
    % -------------------------------
    lefttrials_pc(:,i,:)= currleft_scores;
    righttrials_pc(:,i,:)= currright_scores;
    
    
    % TRY DISTANCE ONLY between PCs averaged over trajectories, rather than
    % tr-by-tr
    % ---------------------------------------
    
    % 1PC only
    left_avgpc1 = meanleft(i,1); % Single value; avg of PC1 for bin i
    right_avgpc1 = meanright(i,1); % If 1 PC, then only 1 value for 1 left and right trial each
    d=dist([left_avgpc1', right_avgpc1']);
    store_avg_DISTPC1(i)=d(1,2); % Will eventually be bin length eg. 11
   
    % 3 PCs
    left_avgpc = meanleft(i,:);
    right_avgpc = meanright(i,:);
    d=dist([left_avgpc', right_avgpc']);  % If 3 PCs, then 3 values for 1 left and right trial each
    store_avg_DISTPC(i)=d(1,2);

    
 
end



%% Do ShufflePCA and get PCA for shuffled trajectory.
% --------------------------------------------------

% d) Distance metric shuffle by taking distance between average shuffle left vs. right, rather
% than individual shuffled trials
store_avg_randDISTPC1=[]; % 1st PC
store_avg_randDISTPC=[]; % All PCs
DistExample=[];
itr=1000;

for shuf = 1:itr    
    % Can generate randperm for each bin separately, or use same randperm
    % for all bins
    
    randtmp = randperm(Ntotaltr);
    randleftidx = randtmp(1:Nlefttr);
    randrightidx = randtmp(Nlefttr+1:end);
        
    for i = 1:Nbins
        X_left = leftvector{i};    % Ntr x Nneu matrix
        X_right = rightvector{i};
        X_combined = [X_left;X_right];
   
        randleftvec = X_combined(randleftidx,:);
        randrightvec = X_combined(randrightidx,:);
        randX_combined = [randleftvec;randrightvec];
        
        [randCOEFF, randSCORES] = pca(randX_combined);
        
        
        % Separate into Random left and right trials, and carry forward
        % Left and Right Shuff PCs
        % randscores is NtrXNpc
        curr_randleftSCORES = randSCORES(1:Nlefttr,:);  %NtrleftXNpc
        curr_randrightSCORES = randSCORES(Nlefttr+1:end,:); %NtrrightXNpc
        currbin_randleftScores_avg = mean(curr_randleftSCORES,1); % mean across trials %1XNpc
        currbin_randrightScores_avg = mean(curr_randrightSCORES,1); % mean across trials %1XNpc
        
       
        % Try getting dist directly from PCs averaged over trials, rather
        % than tr-by-tr below
        % ------------------
        % 1PC
        left_avgpc1=currbin_randleftScores_avg(1); % Single value; avg of PC1 for bin i
        right_avgpc1=currbin_randrightScores_avg(1);
        d=dist([left_avgpc1', right_avgpc1']); % If 1 PC, then only 1 value for 1 left and right trial each
        store_avg_randDISTPC1(shuf,i)=d(1,2); % Will eventually be shufXbin eg. 1000X11

        if shuf==1
            DistExample(i)=d(1,2);
        end
        % 3 PCs
        left_avgpc=currbin_randleftScores_avg(1:3); 
        right_avgpc=currbin_randrightScores_avg(1:3);
        d=dist([left_avgpc', right_avgpc']); % If 3 PCs, then 3 values for 1 left and right trial each
        store_avg_randDISTPC(shuf,i)=d(1,2);
        
               
    end
    
end  % end shuf
   


% IMP
% Get Shuffle DIST PCsvalues Means and CIs - USING AVERAGE ACROSS SHUFFLES
shuffavgDIST_PC1 = mean(store_avg_randDISTPC1,1); % Mean across shuffles, for each timebin  
shuffavgCIDIST_PC1 = [prctile(store_avg_randDISTPC1,5,1);prctile(store_avg_randDISTPC1,95,1)]; 

shuffavgDIST_PC = mean(store_avg_randDISTPC,1); % Mean across shuffles, for each timebin  
shuffavgCIDIST_PC = [prctile(store_avg_randDISTPC,5,1);prctile(store_avg_randDISTPC,95,1)]; 



%% Plot PC distance metric. Real-intertrial distance vs. Shuffled
% -------------------------------------------------------------
oldplot=0;
if oldplot==1
timeaxis = -1000*win(1):binsize_ms:1000*win(2);
timeaxis = timeaxis(1:end-1);


figure(60); hold on; title ('PC Distance Metric; Shuff PCDist; 3PCs black, 1 PC mag');


% This is left vs. right PC distances using average PC - 3 PCs
plot(timeaxis,store_avg_DISTPC,'ko-', 'LineWidth',4,'MarkerSize',16);

% This is shuffle using average left vs right in each shuffle - 3 PCs
plot(timeaxis,shuffavgDIST_PC,'go-', 'LineWidth',4);
plot(timeaxis,shuffavgCIDIST_PC(2,:),'g--', 'LineWidth',2); % 95% CI
plot(timeaxis,shuffavgCIDIST_PC(1,:),'g--', 'LineWidth',2); % 5% CI

% 1 PC - same using average as above
plot(timeaxis,DistExample,'mx-', 'LineWidth',4);
% plot(timeaxis,shuffavgDIST_PC1,'co-', 'LineWidth',4);
% plot(timeaxis,shuffavgCIDIST_PC1(2,:),'c--', 'LineWidth',2); % 95% CI
% plot(timeaxis,shuffavgCIDIST_PC1(1,:),'c--', 'LineWidth',2); % 5% CI


target = find(store_avg_DISTPC>=shuffavgCIDIST_PC(2,:)); % higher than COnfidence interval
Avg_Dis_time = timeaxis(min(target));

% Only look from 100/200ms onward?
% threshbin=4;
% target = find(store_avg_DISTPC(threshbin:end)>=shuffavgCIDIST_PC(2,threshbin:end)) 
% Avg_Dis_time = timeaxis(min(target)+threshbin-1)

legend({'PC distance','Null mean','Null upper 95%CI','Null lower 95%CI','Even-Odd Shuffle'})
end

%% now lets clean this up
timeaxis = -win(1):binsize:win(2);
timeaxis = timeaxis(1:end-1);
regcolors=[242/256 100/256 86/256; 0 162/256 181/256 ];

figure;
title(sprintf('%s Example',region));
% This is left vs. right PC distances using average PC - 3 PCs
hold on;
% This is shuffle using average left vs right in each shuffle - 3 PCs
patch([timeaxis fliplr(timeaxis)], [shuffavgCIDIST_PC(2,:) fliplr(shuffavgCIDIST_PC(1,:))],...
    [.7 .7 .7],'LineStyle','none','FaceAlpha',.5);
plot(timeaxis,shuffavgDIST_PC, 'LineWidth',4,'Color',[.5 .5 .5]);

% 1 PC - same using average as above
plot(timeaxis,DistExample,'k-', 'LineWidth',4);
plot(timeaxis,store_avg_DISTPC, 'LineWidth',4,'Color',regcolors(regionIDX,:));
ylabel('PC Discrimination');
xlabel('Time from odor onset (seconds)');











