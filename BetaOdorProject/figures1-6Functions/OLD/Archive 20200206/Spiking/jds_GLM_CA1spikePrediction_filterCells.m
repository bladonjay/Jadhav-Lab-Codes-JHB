function [gain, s_gain] = jds_GLM_CA1spikePrediction_filterCells(animalprefix,eps)

%%
%Things to incorporate
    %-Need to restrict to cells that pass specific criterion (Active in
    %more than 10 SWR events (Rothschild 2017)
    %Stratification
%%
savedata = 1;

if strcmp(animalprefix,'ER1')
    savedir = ('I:\WT_Singleday\ER1_direct\GLM\');
    ripdir = ('I:\WT_Singleday\ER1_direct\');
    dir = ('I:\WT_Singleday\ER1_direct\');
elseif strcmp(animalprefix,'JS14')
    savedir = ('I:\WT_Singleday\JS14_direct\GLM\');
    ripdir = ('I:\WT_Singleday\JS14_direct\');
    dir = ('I:\WT_Singleday\JS14_direct\');  
elseif strcmp(animalprefix,'KL8')
    savedir = ('I:\WT_Singleday\KL8_direct\GLM\');
    ripdir = ('I:\WT_Singleday\KL8_direct\');
    dir = ('I:\WT_Singleday\KL8_direct\');
elseif strcmp(animalprefix,'JS15')
    savedir = ('I:\WT_Singleday\JS15_direct\GLM\');
    ripdir = ('I:\WT_Singleday\JS15_direct\');
    dir = ('I:\WT_Singleday\JS15_direct\');
elseif strcmp(animalprefix,'JS17')
    savedir = ('I:\WT_Singleday\JS17_direct\GLM\');
    ripdir = ('I:\WT_Singleday\JS17_direct\');
    dir = ('I:\WT_Singleday\JS17_direct\');
elseif strcmp(animalprefix,'JS21')
    savedir = ('I:\WT_Singleday\JS21_direct\GLM\');
    ripdir = ('I:\WT_Singleday\JS21_direct\');
    dir = ('I:\WT_Singleday\JS21_direct\');
end


%%
%-----match neurons across epochs-----%
day = 1; %always single day expt
[ctxidx, hpidx_tmp] = matchidx_acrossep(dir, animalprefix, day); %(tet, cell)
ctxnum = length(ctxidx(:,1));
hpnum_tmp = length(hpidx_tmp(:,1));

%-----create the event matrix during SWRs-----%
spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes

% get ripple time
load(sprintf('%s%srippletime0%d.mat',ripdir,animalprefix,day));

%Get matrix of spikes per cell per SWR event (colums = diff cells, rows =
%diff SWRs) - Same matrix for each PFC cell
gain = [];
s_gain = [];
for ep = eps
    cell2use_idx = zeros(hpnum_tmp,1);
    disp([animalprefix '-epoch-' num2str(ep)])
    
    rip = ripple{day}{ep};
    riptimes(:,1) = rip.starttime;
    riptimes(:,2) = rip.endtime;
    
    rip_starttime = riptimes(:,1)*1000;  % in ms
    
    rip_endtime = riptimes(:,2)*1000;  % in ms
    
    % Find ripples separated by at least 500ms -- sj_getpopulationevents2
    iri = diff(rip_starttime);
    keepidx = [1;find(iri>=500)+1];
    
    riplength = rip_endtime - rip_starttime;
    keepidx2 = find(riplength >= 50);% use the ripple last for more than 50 ms
    keepidx = intersect(keepidx,keepidx2);
    
    riptimes = riptimes(keepidx,:);
    
    %800ms before rip to 600ms before
    preriptimes(:,1) = riptimes(:,1) - 0.2; %start
    preriptimes(:,2) = riptimes(:,1); %end
    
    %start of rip to 200ms after
    postriptimes(:,1) = riptimes(:,1);
    postriptimes(:,2) = riptimes(:,1) + 0.2;
   
    %What peri-SWR ripple to analyze
    hpc_riptimes = postriptimes(:,1);
    hpc_riptimes(:,2) = postriptimes(:,2); %post SWR activity
    
    pfc_riptimes = preriptimes(:,1);
    pfc_riptimes(:,2) = preriptimes(:,2); %pre SWR activity
    
    crit_num = floor(length(hpc_riptimes(:,1))*0.05);
    if crit_num < 10
        crit_num = 10;
    end
    
    if ~isempty(hpc_riptimes)
        for cellcount = 1:hpnum_tmp %get spikes for each cell
            index = [day,ep,hpidx_tmp(cellcount,:)];
            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end
            spikebins = periodAssign(spiketimes, hpc_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
            spk_crit = numel(unique(spikebins)) -    1;
            if spk_crit > crit_num
                goodcell = 1;
            else
                goodcell = 0;
            end
            if goodcell == 1
                cell2use_idx(cellcount,1) = 1;
            else
                cell2use_idx(cellcount,1) = 0;
            end
        end
    end
    
    cell2use_idx = find(cell2use_idx == 1);
    hpidx = hpidx_tmp(cell2use_idx,:);
    hpnum = length(hpidx(:,1));
    
    if ~isempty(pfc_riptimes)
        PFCmatrix = [];
        for cellcount = 1:ctxnum %get spikes for each cell
            index = [day,ep,ctxidx(cellcount,:)] ;
            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end
            %Throw out cells here that are not active for at least 10 SWR
            %events
            spikebins = periodAssign(spiketimes, pfc_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
            if ~isempty(spiketimes)
                validspikes = find(spikebins);
                spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                spikebins = spikebins(validspikes);
            end         
            spikecount = zeros(1,size(pfc_riptimes,1));
            for s = 1:length(spikebins)
                spikecount(spikebins(s)) = spikecount(spikebins(s))+1;
            end
            PFCmatrix = [PFCmatrix spikecount']; %concatenating num spikes per cell, per event
        end
    end
    
    %GET CA1 CELL DATA
    CA1_resp = [];
    if ~isempty(hpc_riptimes)
        for cellcount = 1:hpnum %get spikes for each cell
            index = [day,ep,hpidx(cellcount,:)] ;
            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end
            spikebins = periodAssign(spiketimes, hpc_riptimes(:,[1 2])); %Assign spikes to align with each ripple event (same number = same rip event, number indicates ripple event)
            if ~isempty(spiketimes)
                validspikes = find(spikebins);
                spiketimes = spiketimes(validspikes); %get spike times that happen during ripples
                spikebins = spikebins(validspikes);
            end
            spikecount = zeros(1,size(hpc_riptimes,1));
            for p = 1:length(spikebins)
                spikecount(spikebins(p)) = spikecount(spikebins(p))+1;
            end
            CA1_resp{cellcount}.data = spikecount';
            CA1_resp{cellcount}.CA1cellidx = hpidx(cellcount,:);
        end
    end
    
    %GLM with n-fold cross validation
    mse_CA1 = [];
    for nn = 1:hpnum
        spk_cnt_str = num2str(CA1_resp{nn}.data);
        K = 5;
        cv = cvpartition(spk_cnt_str, 'kfold',K);
        disp(['Cell number ' num2str(nn) ' out of ' num2str(hpnum)])
        mse = zeros(K,1);
        shuf_mse = zeros(K,1);
        for k=1:K
            % training/testing indices for this fold
            trainIdx = cv.training(k);
            testIdx = cv.test(k);
            
            % train GLM model
            CA1mat = CA1_resp{nn}.data;
            CA1mat2 = CA1mat(trainIdx);
            warning('off','all');
            mdl = fitglm(PFCmatrix(trainIdx,:), CA1mat2,'Distribution', 'poisson'); 
            
            % predict regression output
            Y_hat = predict(mdl, PFCmatrix(testIdx,:));
            
            %Do shuffling
            for s = 1:100
                shuf = Y_hat(randperm(length(Y_hat)));
                shuf_err(s) = mean(abs(CA1mat(testIdx) - shuf));
            end
            
            % compute mean squared error
            mse(k) = mean(abs(CA1mat(testIdx) - Y_hat));
            shuf_mse(k) = mean(shuf_err);
            
        end
        mse_CA1{nn}.mse = mse;
        mse_CA1{nn}.mean_mse = mean(mse);
        mse_CA1{nn}.CA1idx = hpidx(nn,:);
        mse_CA1{nn}.mean_shuf_mse = mean(shuf_mse);
        mse_CA1{nn}.shuf_mse = shuf_mse;
        
        %SHUFFLE ORIG DATA AND GET PREDICTION GAIN
        cv2 = cvpartition(spk_cnt_str, 'kfold',K);
        disp('Generating shuffled dataset...')
        for kk=1:K
            disp([num2str(K) ' fold cross validation - fold number ' num2str(kk)])
            trainIdx2 = cv2.training(kk);
            testIdx2 = cv2.test(kk);
            for p = 1:100
                % train GLM model
                shuf_CA1mat = CA1_resp{nn}.data;
                shuf_CA1mat = shuf_CA1mat(randperm(length(shuf_CA1mat)));
                CA1mat2_shuf = shuf_CA1mat(trainIdx2);
                warning('off','all');
                mdl2 = fitglm(PFCmatrix(trainIdx2,:), CA1mat2_shuf,'Distribution', 'poisson'); %Shuffle PFCmat2 to determine shuffled data to get error bars
                
                % predict regression output
                Y_hat_s = predict(mdl2, PFCmatrix(testIdx2,:));
                
                %Do shuffling
                for s = 1:100
                    shuf_shuf = Y_hat_s(randperm(length(Y_hat_s)));
                    shuf_shuf_err(s) = mean(abs(shuf_CA1mat(testIdx2) - shuf_shuf));
                end
                s_gain{ep}{nn}{kk}.data(1,p) = (mean(shuf_shuf_err)/mean(abs(shuf_CA1mat(testIdx2) - Y_hat_s)));
                s_gain{ep}{nn}{kk}.foldnumber = kk;
                s_gain{ep}{nn}{kk}.cellidx = hpidx(nn,:);
                shuf_shuf_err = [];
                if p == 100
                    disp(['Shuffle complete - CA1 cell ' num2str(nn)])
                end
            end
        end
    end
    clear riptimes preriptimes postriptimes hpc_riptimes pfc_riptimes cell2use_idx
  
    pred_gain = [];
    
    for mm = 1:length(mse_CA1)
        mse = mean(mse_CA1{mm}.mse);
        shuf_mse = mean(mse_CA1{mm}.shuf_mse);
        pg = (shuf_mse/mse);
        pred_gain = [pred_gain pg];
    end
    
    m_pred_gain(ep/2) = mean(pred_gain);
%     m_pred_gain(ep) = mean(pred_gain);
    gain{ep}.meanGain = mean(pred_gain);
    gain{ep}.allCellGain = pred_gain;
    gain{ep}.animal = animalprefix;
    gain{ep}.epoch = ep;
end

if savedata == 1
   save(sprintf('%s%s_SWR_GLMfiltStrat_PFCpredCA1_PRE200-0toPOST0-200_changeErrorCalc_PG%02d.mat', savedir,animalprefix,day), 'gain');
   save(sprintf('%s%s_SWR_GLMfiltStrat_PFCpredCA1_PRE200-0toPOST0-200_changeErrorCalc_shuffled_PG%02d.mat', savedir,animalprefix,day), 's_gain');
end

disp([animalprefix ' processing complete'])

