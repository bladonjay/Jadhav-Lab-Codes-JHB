%cs_GLMtrialOutcomePrediction

%predict left vs right, or correct vs incorrect based on spiking during
%sniffing.

% Do separately for each day? How do we combine across days/animals (different
% cells recorded)

function [fract_correct, fract_correct_shuff] = cs_predictChoice_GLM(animal, day, region, window,celltype)
topDir = cs_setPaths();
animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
fullwin = 1;

%get matrix of spikes per cell per trial (colums = diff cells, rows =
%diff trials)
%assign trial IDs to each trial
%run glm - train on subset of data, binomial distribution
%validation?

%%

%-----create the event matrix-----%
spikes = loaddatastruct(animDir, animal, 'spikes', day); % get spikes
cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
switch celltype
    case 'selective'
        filt = ['strcmp($area,''',region,''') && $numspikes>0 && strcmp($type,''pyr'') && (strcmp($selectivity,''rightSelective'') || strcmp($selectivity,''leftSelective''))'];
    case 'all'
        filt = ['strcmp($area,''',region,''') && $numspikes>0 && strcmp($type,''pyr'')'];
end
cells = evaluatefilter(cellinfo{day},filt);
cells = unique(cells(:,[2 3]),'rows');
numcells = size(cells,1);

% get get trial times and ids
odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers', day);

eps = cs_getRunEpochs(animDir, animal, 'odorplace',day);
eps = eps(:,2);
%Get matrix of spikes per cell per trial (colums = diff cells, rows =
%diff trials)
gain = [];
s_gain = [];

Cellmatrix = [];
Cellmatrix_fullwin = [];
for c = 1:numcells %get spikes for each cell
    goodeps = cs_findGoodEpochs(spikes{day},{'data'},cells(c,:)); %dont include cells with no spikes during run eps
    spikecount_day = [];
    spikecount_fullwin = [];
    if ~isempty(intersect(goodeps, eps)) %cell should spike during at least one run epoch
        for ep = eps'
            trials = [];
            index = [day,ep,cells(c,:)] ;
            trials(:,1) = odorTriggers{day}{ep}.allTriggers;
            trials(:,2) = trials(:,1) + window;

            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end
            
            %Assign spikes to align with each trial (same number = same trial, number indicates trial)
            spiketrials = periodAssign(spiketimes, trials(:,[1 2])); 
            
            spiketrials = nonzeros(spiketrials);
            
            spikecount = zeros(1,size(trials,1));
            for s = 1:length(spiketrials)
                spikecount(spiketrials(s)) = spikecount(spiketrials(s))+1;
            end
            spikecount_day = [spikecount_day; spikecount']; %concatenating num spikes per cell, per event
            
            
            %get the spikecount in the full time window for odor sampling.
            %Use this to determine whether cell should be excluded based on
            %low FR. Otherwise, don't have the same number of cells across
            %time windows
            
            fulltrials = [trials(:,1), trials(:,1) + fullwin];
            
            %Assign spikes to align with each trial (same number = same trial, number indicates trial)
            spiketrials_full = periodAssign(spiketimes, fulltrials(:,[1 2])); 
            
            spiketrials_full = nonzeros(spiketrials_full);
            
            spikecount_full = zeros(1,size(fulltrials,1));
            for s = 1:length(spiketrials_full)
                spikecount_full(spiketrials_full(s)) = spikecount_full(spiketrials_full(s))+1;
            end
            spikecount_fullwin = [spikecount_fullwin; spikecount_full']; %concatenating num spikes per cell, per event
            
        end
        Cellmatrix_fullwin = [Cellmatrix_fullwin spikecount_fullwin]; %sum spikes across epochs, store in a row of array for each cell
        Cellmatrix = [Cellmatrix spikecount_day]; %concatenating cells
    end
end

 %Throw out cells here that are not active on at least 10 trials during
 %full NP window
Cellmatrix_fullwin(Cellmatrix_fullwin ~= 0) = 1;
activetrials = sum(Cellmatrix_fullwin,1);
validcells = activetrials>= 10;
Cellmatrix = Cellmatrix(:,validcells);


%If fewer than 3 cells, don't try GLM 
if size(Cellmatrix,2) < 3
    fract_correct = [];
    fract_correct_shuff = [];
    return
end
%GET TRIAL DATA
trialIDs = [];
for ep = eps'
    [cl, cr, il, ir] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
    trials = zeros(length(odorTriggers{day}{ep}.allTriggers),1);
    
    %set trials where animal ultimately chose left reward well to 1, regardless
    %of correct/incorrect
    trials([cl;ir]) = 1;
    trialIDs = [trialIDs;trials];
end

numtrials = length(trialIDs);

fract_correct = [];

K = 5;
cv = cvpartition(numtrials, 'kfold',K);
%        
for k=1:K
    % training/testing indices for this fold
    trainIdx = cv.training(k);
    testIdx = cv.test(k);
    
    % train GLM model
    
    warning('off','all');
    %disp('Fitting GLM')
    mdl = fitglm(Cellmatrix(trainIdx,:), trialIDs(trainIdx),'Distribution', 'binomial');
    
    % predict regression output
    Y_hat = predict(mdl, Cellmatrix(testIdx,:));
    
pred = round(Y_hat); %glm outputs very small values instead of zeros (binomial) for some reason, round to zero
ids = trialIDs(testIdx);
fract_correct(k) = sum(pred == ids)/length(pred);

end
fract_correct = mean(fract_correct);

%Do shuffling
for s = 1:100
    if mod(s,10) == 0 || s == 1
    %disp(['Doing shuffle- iteration ', num2str(s)])
    end
    
    shuff_trialIDs = trialIDs(randperm(length(trialIDs)));
    mdl = fitglm(Cellmatrix(trainIdx,:), shuff_trialIDs(trainIdx),'Distribution', 'binomial');
    Y_hat_shuff = predict(mdl, Cellmatrix(testIdx,:));
    pred_shuff = round(Y_hat_shuff);
    ids_shuff = shuff_trialIDs(testIdx);
    
    fract_correct_shuff(s) = sum(pred_shuff == ids_shuff)/length(pred_shuff);
end
fract_correct_shuff = mean(fract_correct_shuff);

