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

region = 'CA1'; 
%animal = {'CS34'}; 
animal = 'CS33'; animno = 2; day = 1;
dataDir = '/Users/Shantanu/data25/OLF_CS/Data/';
figdir = '/Users/Shantanu/data25/OLF_CS/Data/sj_Figures';

animalDir = [dataDir,animal,'_direct/'];
load([animalDir,animal,'cellinfo.mat'])
daystr = getTwoDigitNumber(day);
load([animalDir,animal,'nosepokeWindow',daystr,'.mat'])

win = [0.2 0.9]; % Implies from -0.2 to 1
binsize=0.1;
binsize_ms=1000*binsize;
selectiveonly=0;
timeaxis = -1000*win(1):binsize_ms:1000*win(2);

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
allbins = [-0.5:binsize:1.5-binsize];
goodbins = [win1:binsize:win(2)-binsize]; 
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
SCORE = SCORES;
% COEFF returns N_neurons prncipal components in columns, with each PC
% having Nneurons coefficients along original dimensions

% SCORES are the projections: Original data in Ntrials projected along
% the new N_neu dimensions - we will look at the first 3




%%  Shuffle PCA firing rate
Nneu = size(lefttrials,1)
Nlefttr = length(leftidx)
Nrighttr = length(rightidx)
Ntotaltr = Nlefttr + Nrighttr;
itr=1000;



for shuf = 1:itr
    randtmp = randperm(Ntotaltr);
    randleftidx = randtmp(1:Nlefttr);
    randrightidx = randtmp(Nlefttr+1:end);
    
    randleftfr = LRcombined(randleftidx,:);
    randrightfr = LRcombined(randrightidx,:);
    randLRcombined = [randleftfr;randrightfr];
    
    [randCOEFF, randSCORES] = pca(randLRcombined);
    randfrSCORES_PC1(shuf,:) = randSCORES(:,1);
    randfrSCORES_PC2(shuf,:) = randSCORES(:,2);
    randfrSCORES_PC3(shuf,:) = randSCORES(:,3);
    randfrleftSCORES_PC1(shuf,:) = randSCORES(1:Nlefttr,1);
    randfrleftSCORES_PC2(shuf,:) = randSCORES(1:Nlefttr,2);
    randfrleftSCORES_PC3(shuf,:) = randSCORES(1:Nlefttr,3);
    randfrrightSCORES_PC1(shuf,:) = randSCORES(Nlefttr+1:end,1);
    randfrrightSCORES_PC2(shuf,:) = randSCORES(Nlefttr+1:end,2);
    randfrrightSCORES_PC3(shuf,:) = randSCORES(Nlefttr+1:end,3);
end
    



%% Firing rate Figures

% figure('color','w'); title (['PC Firing Rates']);
% scatter3(SCORE(leftidx,1),SCORE(leftidx,2),SCORE(leftidx,3),100,...
%                     'filled','MarkerFaceColor','r',...
%                         'MarkerEdgeColor','w');                    
% hold on,
% scatter3(SCORE(rightidx,1),SCORE(rightidx,2),SCORE(rightidx,3),100,...
%                     'filled','MarkerFaceColor','b',...
%                         'MarkerEdgeColor','w');
% 
% % Plot Shuffle Means
% plot3(mean(randfrleftSCORES_PC1,1),mean(randfrleftSCORES_PC2,1),mean(randfrleftSCORES_PC3,1),'mx');
% plot3(mean(randfrrightSCORES_PC1,1),mean(randfrrightSCORES_PC2,1),mean(randfrrightSCORES_PC3,1),'bx');                   
%                     
%                     
%                     
%                      
% figure;
% subplot(2,1,1);
% bar(EXPLAINED);
% ylabel('Percentage variance explained');
% xlabel('Dimension number');
% box off;
% axis([0 3.5 0 100]);
% 
% subplot(2,1,2);
% bar(cumsum(EXPLAINED));
% ylabel('Percentage total variance explained using components 1-N');
% xlabel('Dimension number');
% box off;
% axis([0 3.5 0 100]);
%                     





%% 

% PCA on firing rate vectors
% ------------------------

leftvector=[]; rightvector=[]; %For each timebin, save NtrialXNnenuron firing matrix

Nneu = size(lefttrials,1)
Nlefttr = length(leftidx)
Nrighttr = length(rightidx)
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


% 2022 - The following is not used anymore

% % For storing Distance Metric in PC values
% % ---------------------------------------
% store_left_DISTPC=[]; % Intra-left trial distance in each bin
% store_right_DISTPC=[]; % Intra-right trial distance in each bin
% store_leftmean_DISTPC=[]; % mean intra-left DIST in each bin
% store_rightmean_DISTPC=[]; % mean intra-right DIST in each bin
% store_leftstd_DISTPC=[]; % std intra-left DIST in each bin
% store_rightstd_DISTPC=[]; % std intra-right DIST in each bin
% store_leftCI_DISTPC=[]; % 5% and 95% CI intra-left DIST in each bin
% store_rightCI_DISTPC=[]; % 5% and 95% CIstd intra-right DIST in each bin
% 
% store_DISTPC = []; % Inter left-right trial distance in each bin
% store_mean_DISTPC=[]; % mean Inter left-right DIST in each bin
% store_std_DISTPC=[]; % std Inter left-right DIST in each bin
% store_CI_DISTPC=[]; % 5% and 95% CI Inter left-right DIST in each bin



for i = 1:Nbins
    X_left = leftvector{i};
    X_right = rightvector{i};
    
    LRcombined = [X_left; X_right];
    [COEFF, SCORES, LATENT, TSQUARED, EXPLAINED, MU] = pca(LRcombined);
    SCORE = SCORES;
    
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
    
    
    
    % % 2022: NOTE - The fOllowing is NOT USED anymore
    
%     % Get Distance Metric between 1st 3 PCs
%     % ------------------------------------
%     % All trials
%     X_comb = SCORES(:,1:3); % Ntrtotal x 3 PCs
%     X_dist = dist(X_comb');  % Xcomb' has Ntrleft+Ntrright columns
%     X_dist_inter = X_dist(1:Nlefttr,Nlefttr+1:end);
%     X_dist_left = X_dist(1:Nlefttr,1:Nlefttr);
%     X_dist_right = X_dist(Nlefttr+1:end,Nlefttr+1:end);
%     
%      % Left trials
%     % ------------
%     X_dist_left = triu(X_dist_left,1);
%     %dist_currbin_leftvec = X_dist_left(X_dist_left~=0);
%     dist_currbin_leftvec = X_dist_left(~tril(ones(size(X_dist_left))));
%     store_left_DISTPC{i} = dist_currbin_leftvec;
%     store_leftmean_DISTPC(i) = mean(dist_currbin_leftvec);
%     store_leftstd_DISTPC(i) = std(dist_currbin_leftvec);
%     store_leftCI_DISTPC(i,:) = [prctile(dist_currbin_leftvec,5),prctile(dist_currbin_leftvec,95)];
%     
%     % Right trials
%     % ------------
%     X_dist_right = triu(X_dist_right,1);
%     %dist_currbin_rightvec = X_dist_right(X_dist_right~=0);
%     dist_currbin_rightvec = X_dist_right(~tril(ones(size(X_dist_right))));
%     store_right_DISTPC{i} = dist_currbin_rightvec;
%     store_rightmean_DISTPC(i) = mean(dist_currbin_rightvec);
%     store_rightstd_DISTPC(i) = std(dist_currbin_rightvec);
%     store_rightCI_DISTPC(i,:) = [prctile(dist_currbin_rightvec,5),prctile(dist_currbin_rightvec,95)];
%     
%     % Inter-trial Distance, Left-Right
%     % --------------------------------
%    
%     X_dist_inter = triu(X_dist_inter,1);
%     %dist_currbin_vec = X_dist_inter(X_dist_inter~=0);
%     dist_currbin_vec = X_dist_inter(~tril(ones(size(X_dist_inter))));
%     store_DISTPC{i} = dist_currbin_vec;
%     store_mean_DISTPC(i) = mean(dist_currbin_vec);
%     store_std_DISTPC(i) = std(dist_currbin_vec);
%     store_CI_DISTPC(i,:) = [prctile(dist_currbin_vec,5),prctile(dist_currbin_vec,95)];
    
    
 
end

% figure; hold on; title (['Explained Variance in each time bin']);
% bar(store_Explained3);
% ylabel('Percentage variance explained by 3 PCs');
% xlabel('Time Bin number');
% box off;
% axis([0 13.5 0 100]);



% % PLOT DISTANCE METRIC using 1st 3 PCs
% % ------------------------------------------
% 
 timeaxis = -1000*win(1):binsize_ms:1000*win(2);
 timeaxis = timeaxis(1:end-1);

% % 2022: NOTE - The fOllowing is NOT USED anymore
% figure(6); hold on; title ('PC Distance Metric');
% plot(timeaxis,store_leftmean_DISTPC,'ro-', 'LineWidth',4);
% plot(timeaxis,store_leftmean_DISTPC+store_leftstd_DISTPC,'r--', 'LineWidth',2);
% plot(timeaxis,store_leftmean_DISTPC-store_leftstd_DISTPC,'r--', 'LineWidth',2);
% 
% plot(timeaxis,store_rightmean_DISTPC,'bx-', 'LineWidth',4);
% plot(timeaxis,store_rightmean_DISTPC+store_rightstd_DISTPC,'b--', 'LineWidth',2);
% plot(timeaxis,store_rightmean_DISTPC-store_rightstd_DISTPC,'b--', 'LineWidth',2);
% 
% 
% plot(timeaxis,store_mean_DISTPC,'yx-', 'LineWidth',2);
% plot(timeaxis,store_mean_DISTPC+store_std_DISTPC,'y--', 'LineWidth',2);
% plot(timeaxis,store_mean_DISTPC-store_std_DISTPC,'y--', 'LineWidth',2);
% 
% 
% plot(timeaxis,store_avg_DISTPC,'kx-', 'LineWidth',4);








%% Do ShufflePCA and get PCA for shuffled trajectory.
% ----------------------------------------------------

% 2022: NOT USES Options a,b,c, BELOW ANYMORE. Only d


% a) Don't divide into left and right
% b) Divide into left and right and get 2 shuffled Left/Right trajectories
% - can compute distance between them. Better done in real data.

% %a) Total trials - don't divide
% store_randScores=[]; % shufnum X binnum X Ntr X NneuNpc
% store_randScores_PC1=[]; % shufnum X binnum X Ntr
% store_randScores_PC2=[]; % shufnum X binnum X Ntr
% store_randScores_PC3=[]; % shufnum X binnum X Ntr
% 
% store_randScores_avgtr=[]; %Shuf X binnum X Nneu/Npc - avg across trials.
% store_randScores_avgtr_PC1=[];
% store_randScores_avgtr_PC2=[];
% store_randScores_avgtr_PC3=[];
% 
% %b) Divide into left and right
% storeleft_randScores_avgtr=[]; %Shuf X binnum X Nneu/Npc - avg across trials.
% storeleft_randScores_avgtr_PC1=[];
% storeleft_randScores_avgtr_PC2=[];
% storeleft_randScores_avgtr_PC3=[];
% 
% storeright_randScores_avgtr=[]; %Shuf X binnum X Nneu/Npc - avg across trials.
% storeright_randScores_avgtr_PC1=[];
% storeright_randScores_avgtr_PC2=[];
% storeright_randScores_avgtr_PC3=[];


% c) Get distance metric shuffle for 1st 3 PCs
store_randDISTPC = []; % SHUF Inter left-right trial distance in each bin
store_mean_randDISTPC=[]; % mean SHUF Inter left-right DIST in each bin for 1st 3 PCs
%store_std_randDISTPC=[]; % std Inter left-right DIST in each bin
%store_CI_randDISTPC=[]; % 5% and 95% CI Inter left-right DIST in each bin


% d) Distance metric shuffle by taking distance between average shuffle left vs. right, rather
% than individual shuffled trials
store_avg_randDISTPC1=[]; % 1st PC
store_avg_randDISTPC=[]; % All PCs


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
        
        
        % REMOVE FOR AVG DIS
        
        % randscores is NtrXNpc
%         %  REMOVE FOR AVG DIS
%         store_randScores(shuf,i,:,:)=randSCORES;  %shufXbinXNtrXNneu
%         store_randScores_PC1(shuf,i,:)=randSCORES(:,1);  %shufXbinXNtr
%         store_randScores_PC2(shuf,i,:)=randSCORES(:,2);  %shufXbinXNtr
%         store_randScores_PC3(shuf,i,:)=randSCORES(:,3);  %shufXbinXNtr
        
        
        % Don't need to carry forward trial numbers in shuffle. So get
        % average scores across trials
        
%         % REMOVE FOR AVG DIS
%         currbin_randScores_avg = mean(randSCORES,1); % Mean across trials. Should be vector with length(Nneu/NPc)
%         store_randScores_avgtr(shuf,i,:) = currbin_randScores_avg; %shufXbinXNneu/Npc
%         store_randScores_avgtr_PC1(shuf,i,:) = currbin_randScores_avg(1); % shufXbin  eg. 1000X11
%         store_randScores_avgtr_PC2(shuf,i,:) = currbin_randScores_avg(2); % shufXbin
%         store_randScores_avgtr_PC3(shuf,i,:) = currbin_randScores_avg(3); % shufXbin
        
        
        
        % Separate into Random left and right trials, and carry forward
        % Left and Right Shuff PCs
        % randscores is NtrXNpc
        curr_randleftSCORES = randSCORES(1:Nlefttr,:);  %NtrleftXNpc
        curr_randrightSCORES = randSCORES(Nlefttr+1:end,:); %NtrrightXNpc
        currbin_randleftScores_avg = mean(curr_randleftSCORES,1); % mean across trials %1XNpc
        currbin_randrightScores_avg = mean(curr_randrightSCORES,1); % mean across trials %1XNpc
        
        
%         % REMOVE FOR AVG DIS
%         storeleft_randScores_avgtr(shuf,i,:) = currbin_randleftScores_avg; %shufXbinXNpc
%         storeleft_randScores_avgtr_PC1(shuf,i,:) = currbin_randleftScores_avg(1); % shufXbin  eg. 1000X11
%         storeleft_randScores_avgtr_PC2(shuf,i,:) = currbin_randleftScores_avg(2); % shufXbin
%         storeleft_randScores_avgtr_PC3(shuf,i,:) = currbin_randleftScores_avg(3); % shufXbin
%         
%         storeright_randScores_avgtr(shuf,i,:) = currbin_randrightScores_avg; %shufXbinXNpc
%         storeright_randScores_avgtr_PC1(shuf,i,:) = currbin_randrightScores_avg(1); % shufXbin  eg. 1000X11
%         storeright_randScores_avgtr_PC2(shuf,i,:) = currbin_randrightScores_avg(2); % shufXbin
%         storeright_randScores_avgtr_PC3(shuf,i,:) = currbin_randrightScores_avg(3); % shufXbin
        
        
        
        
        
        
        % Try getting dist directly from PCs averaged over trials, rather
        % than tr-by-tr below
        % ------------------
        % 1PC
        left_avgpc1=currbin_randleftScores_avg(1); % Single value; avg of PC1 for bin i
        right_avgpc1=currbin_randrightScores_avg(1);
        d=dist([left_avgpc1', right_avgpc1']); % If 1 PC, then only 1 value for 1 left and right trial each
        store_avg_randDISTPC1(shuf,i)=d(1,2); % Will eventually be shufXbin eg. 1000X11

        % 3 PCs
        left_avgpc=currbin_randleftScores_avg(1:3); 
        right_avgpc=currbin_randrightScores_avg(1:3);
        d=dist([left_avgpc', right_avgpc']); % If 3 PCs, then 3 values for 1 left and right trial each
        store_avg_randDISTPC(shuf,i)=d(1,2);
        
        
        
%         % REMOVE FOR AVG DIS
% 
%         % Shuffled PC distance metric
%         % Get Distance Metric between ist 3 PCs or 1 PC
%         % ------------------------------------
%         % All trials
%         X_comb = randSCORES(:,1:3); % Ntrtotal x 3 PCs
%         %X_comb = randSCORES(:,1); % Ntrtotal x 1 PC
%         X_dist = dist(X_comb');  % Xcomb' has Ntrleft+Ntrright columns
%         X_dist_inter = X_dist(1:Nlefttr,Nlefttr+1:end);
%         
%         
%         % Random Inter-trial Distance, Left-Right
%         % --------------------------------
%         
%         X_dist_inter = triu(X_dist_inter,1);
%         %dist_currbin_vec = X_dist_inter(X_dist_inter~=0);
%         dist_currbin_vec = X_dist_inter(~tril(ones(size(X_dist_inter))));
%         store_randDISTPC(shuf,i,:) = dist_currbin_vec;
%         store_mean_randDISTPC(shuf,i) = mean(dist_currbin_vec);
%         %store_std_randDISTPC(i) = std(dist_currbin_vec);
%         %store_CI_randDISTPC(i,:) = [prctile(dist_currbin_vec,5),prctile(dist_currbin_vec,95)];
        
        % Getting these distances one-by-one through matrix, rather than
        % using dist on whole matrix. For 1 PC
%         for a=1:43; Xdistinter(a) = mean(X_dist(a,44:end)); end
%         mean(Xdistinter);
%         for a=1:43; Xdistleft(a) = mean(X_dist(a,1:43)); end
%         mean(Xdistleft)
%         for a=44:84; Xdistright(a-43) = mean(X_dist(a,44:84)); end
%         mean(Xdistright)
    end
    
end  % end shuf
   



% %         % REMOVE FOR AVG DIS
% 
% % Get Shuffle values Means and CIs
% shuffmean_PC1 = mean(store_randScores_avgtr_PC1,1); % Mean shuf value across timebins
% %shuffCI_PC1 = [shuffmean_PC1 - prctile(store_randScores_avgtr_PC1,5,1); shuffmean_PC1 + prctile(store_randScores_avgtr_PC1,95)];
% shuffCI_PC1 = [prctile(store_randScores_avgtr_PC1,5,1); prctile(store_randScores_avgtr_PC1,95,1)];
% 
% shuffmeanleft_PC1 = mean(storeleft_randScores_avgtr_PC1,1); % Mean shuf value across timebins
% %shuffCIleft_PC1 = [shuffmeanleft_PC1 - prctile(storeleft_randScores_avgtr_PC1,5,1); shuffmeanleft_PC1 + prctile(storeleft_randScores_avgtr_PC1,95)];
% shuffCIleft_PC1 = [prctile(storeleft_randScores_avgtr_PC1,5,1); prctile(storeleft_randScores_avgtr_PC1,95,1)];
% 
% shuffmeanright_PC1 = mean(storeright_randScores_avgtr_PC1,1); % Mean shuf value across timebins
% %shuffCIright_PC1 = [shuffmeanright_PC1 - prctile(storeright_randScores_avgtr_PC1,5,1); shuffmeanright_PC1 + prctile(storeright_randScores_avgtr_PC1,95)];
% shuffCIright_PC1 = [prctile(storeright_randScores_avgtr_PC1,5,1); prctile(storeright_randScores_avgtr_PC1,95,1)];
% 
% % Get Shuffle DIST PCsvalues Means and CIs USING ALL SHUFFLED TRIALS
% shuffmeanDIST_PC = mean(store_mean_randDISTPC,1); % Mean across shuffles, for each timebin  
% shuffCIDIST_PC = [prctile(store_mean_randDISTPC,5,1);prctile(store_mean_randDISTPC,95,1)]; 


% IMP
% Get Shuffle DIST PCsvalues Means and CIs - USING AVERAGE ACROSS SHUFFLES
shuffavgDIST_PC1 = mean(store_avg_randDISTPC1,1); % Mean across shuffles, for each timebin  
shuffavgCIDIST_PC1 = [prctile(store_avg_randDISTPC1,5,1);prctile(store_avg_randDISTPC1,95,1)]; 

shuffavgDIST_PC = mean(store_avg_randDISTPC,1); % Mean across shuffles, for each timebin  
shuffavgCIDIST_PC = [prctile(store_avg_randDISTPC,5,1);prctile(store_avg_randDISTPC,95,1)]; 




%% Plot PC distance metric. Real-intertrial distance vs. Shuffled
% -------------------------------------------------------------


figure(60); hold on; title ('PC Distance Metric; Shuff PCDist; 3PCs black, 1 PC mag');

% This is left vs. right PC distances using all tr-to-tr distances
%plot(timeaxis,store_mean_DISTPC,'kx-', 'LineWidth',2);
%plot(timeaxis,store_mean_DISTPC+store_std_DISTPC,'k--', 'LineWidth',2);
%plot(timeaxis,store_mean_DISTPC-store_std_DISTPC,'k--', 'LineWidth',2);

% This is shuffle using all tr-to-tr distances
%plot(timeaxis,shuffmeanDIST_PC,'go-', 'LineWidth',2);
%plot(timeaxis,shuffCIDIST_PC(2,:),'g--', 'LineWidth',2); % 95% CI
%plot(timeaxis,shuffCIDIST_PC(1,:),'g--', 'LineWidth',2); % 5% CI

% This is left vs. right PC distances using average PC - 3 PCs
plot(timeaxis,store_avg_DISTPC,'ko-', 'LineWidth',4,'MarkerSize',16);

% This is shuffle using average left vs right in each shuffle - 3 PCs
plot(timeaxis,shuffavgDIST_PC,'go-', 'LineWidth',4);
plot(timeaxis,shuffavgCIDIST_PC(2,:),'g--', 'LineWidth',2); % 95% CI
plot(timeaxis,shuffavgCIDIST_PC(1,:),'g--', 'LineWidth',2); % 5% CI

% 1 PC - same using average as above
% plot(timeaxis,store_avg_DISTPC1,'mx-', 'LineWidth',4);
% plot(timeaxis,shuffavgDIST_PC1,'co-', 'LineWidth',4);
% plot(timeaxis,shuffavgCIDIST_PC1(2,:),'c--', 'LineWidth',2); % 95% CI
% plot(timeaxis,shuffavgCIDIST_PC1(1,:),'c--', 'LineWidth',2); % 5% CI


target = find(store_avg_DISTPC>=shuffavgCIDIST_PC(2,:)); % higher than COnfidence interval
Avg_Dis_time = timeaxis(min(target));

% Only look from 100/200ms onward?
% threshbin=4;
% target = find(store_avg_DISTPC(threshbin:end)>=shuffavgCIDIST_PC(2,threshbin:end)) 
% Avg_Dis_time = timeaxis(min(target)+threshbin-1)

Nneu

keyboard;















%%         % REMOVE FOR AVG DIS



% % %  PLOTS FOR 1-d PCA
% % % -----------------
% % 
% % % 1d plot, with sem/shuffle/95% interval
% % % -----------------------
% PCnum=1;
% % 
% % % PLOT mean and sem
% % % -----------------
% figure(10); hold on; title ('1-d PC with sem; Shuff PCs');
% plot(timeaxis,meanleft(:,PCnum),'ro-', 'LineWidth',4);
% plot(timeaxis,meanleft(:,PCnum)+errleft(:,PCnum),'r--', 'LineWidth',2);
% plot(timeaxis,meanleft(:,PCnum)-errleft(:,PCnum),'r--', 'LineWidth',2);
% 
% plot(timeaxis,meanright(:,PCnum),'bx-', 'LineWidth',4);
% plot(timeaxis,meanright(:,PCnum)+errright(:,PCnum),'b--', 'LineWidth',2);
% plot(timeaxis,meanright(:,PCnum)-errright(:,PCnum),'b--', 'LineWidth',2);
% % 
% 
% %% 3d plot
% figure(99); hold on;
% plot3(meanleft(:,1),meanleft(:,2),meanleft(:,3),'ro-', 'LineWidth',4);
% plot3(meanright(:,1),meanright(:,2),meanright(:,3),'bx-', 'LineWidth',4);
% 
% PCnum
% 
% keyboard;





% % Single trials
% % -----------------------

% % Pick 10 random trials
% % randtmp = randperm(Nlefttr);
% % %Nlefttr=Nlefttr(randtmp(1:10));
% % %for tr = 1:Nlefttr
% % for tr = 1:3
% %     currvec = squeeze(lefttrials_pc(randtmp(tr),:,:));
% %     plot3(currvec(:,1),currvec(:,2),currvec(:,3),'r-');
% % end
% % 

% % randtmp = randperm(Nrighttr);
% % %Nrighttr=Nrighttr(randtmp(1:10));
% % %for tr = 1:Nrighttr
% % for tr = 1:3
% %     currvec = squeeze(righttrials_pc(randtmp(tr),:,:));
% %     plot3(currvec(:,1),currvec(:,2),currvec(:,3),'b-');
% % end
% 

% % 2d plot
% figure(100); hold on;
% plot(meanleft(:,1),meanleft(:,2),'ro-', 'LineWidth',4);
% plot(meanright(:,1),meanright(:,2),'bx-', 'LineWidth',4);
% 
% % Single trials - all
% % %for tr = 1:Nlefttr
% % for tr = 1:3
% %     currvec = squeeze(lefttrials_pc(tr,:,:));
% %     plot(currvec(:,1),currvec(:,2),'r-');
% % end
% % 
% % %for tr = 1:Nrighttr
% % for tr = 1:3
% %     currvec = squeeze(righttrials_pc(tr,:,:));
% %     plot(currvec(:,1),currvec(:,2),'b-');
% % end
% 
% %% 
%  
% %% Single trials - pick randomly
% randtmp = randperm(Nlefttr);
% %Nlefttr=Nlefttr(randtmp(1:10));
% %for tr = 1:Nlefttr
% for tr = 1:3
%     currvec = squeeze(lefttrials_pc(randtmp(tr),:,:));
%     plot(timeaxis,currvec(:,1),'m.-','LineWidth',1);
% end
% 
% randtmp = randperm(Nrighttr);
% %Nrighttr=Nrighttr(randtmp(1:10));
% %for tr = 1:Nrighttr
% for tr = 1:3
%     currvec = squeeze(righttrials_pc(randtmp(tr),:,:));
%     plot(timeaxis,currvec(:,1),'c.-','LineWidth',1);
% end
% 
keyboard;



