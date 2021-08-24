function spatial_reconstruction_CA1PFC_singleday(animalprefix,day,ep,cellcountthresh,wellcutoff,savedata,savedir)
%---------------------------------------------------------------%
%  This is the function for replay detection using CA1 spikes   %
%  during SWRs                                                  %
%  -- Wenbo Tang (Sep 13, 2019)                                 %
%---------------------------------------------------------------%

% INPUTS:
%
%    animalprefix = animal prefix.
%    day = experimental day.
%    ep = epoch.
%    cellcountthresh = mininum number of cells active for considering as a
%                      cadidate event, usually = 5
%    wellcutoff = region excluded, cm from start and end
%    savedata = save results, 1 = save, 0 = not save
%    savedir = directory for saving results

%%
%---- add the codes to paths ---%
addpath(genpath('/Users/wenbotang/Src_Matlab'))
%---- set parameters ----%
tBinSz = 10; %default temporal bin in ms used for replay detection, hard coded

nstd=round(20/tBinSz); % 20ms gaussian kernal for reactivation smooth
g1 = gaussian(nstd, 5*nstd+1);
%%
% set animal directory
if strcmp(animalprefix,'ER1')
    dir = 'D:\SingledayExp\ER1_NEW_direct2\';
    exclude_list = [1 1;1 2]; % exclude interneurons
elseif strcmp(animalprefix,'KL8')
    dir = 'D:\SingledayExp\KL8_direct\';
    exclude_list = [];
elseif strcmp(animalprefix,'JS14')
    dir = 'D:\SingledayExp\JS14_direct\';
    exclude_list = [6,1;8,3;23,1];
elseif strcmp(animalprefix,'JS15')
    dir = 'D:\SingledayExp\JS15_direct\';
    exclude_list = [5,4;8,6];
elseif strcmp(animalprefix,'JS17')
    dir = 'D:\SingledayExp\JS17_direct\';
    exclude_list = [6,2;6,4;6,6;7,5;10,1;11,2;23,2;23,3;23,4];
elseif strcmp(animalprefix,'JS21')
    dir = 'D:\SingledayExp\JS21_direct\';
    exclude_list = [6,5;21,2;25,2;25,3];
end
%%
% load previous file to add new result
if ep > 2
    load(sprintf('%s%sreactivationtraj_hpctx_%02d.mat', savedir,animalprefix,day));
end
%%
%-----match neurons across epochs-----%
[ctxidx, hpidx] = matchidx_acrossep_singleday(dir, animalprefix, day,exclude_list); %(tet, cell)

ctxnum = length(ctxidx(:,1));
hpnum = length(hpidx(:,1));
%%
%-----create the ratemaps [nPosBin x nHPCells]-----%
rm = []; % ratemap matrix
pm = []; % position matrix
tm = []; % track matrix
cellidxm_hp = []; % CA1 cell index
cellidxm_ctx = [];% PFC cell index

load(sprintf('%s%slinfields0%d.mat',dir,animalprefix,day)); % get linearized place fields
% CA1 cell loop
for i = 1:hpnum
      cind = hpidx(i,:);
      if (length(linfields{day}{eprun})>= cind(1))
            if  (length(linfields{day}{eprun}{cind(1)})>= cind(2))
                linfield1 = linfields{day}{eprun}{cind(1)}{cind(2)};
            else 
                linfield1 =[];
            end
      else
            linfield1=[];
      end
      
      if ~isempty(linfield1)
           linfield_hp = [];
           lintrack_hp = [];
           pos_hp = [];
           % 4 different trajectory types
           for track = 1:4
                temp1 = linfield1{track};
                pos1 = temp1(:,1);
                lintrack1 = ones(size(pos1))*track;
                occnormrate1 = temp1(:,5);
                linfield_hp = [linfield_hp;occnormrate1];
                pos_hp = [pos_hp;pos1];
                lintrack_hp = [lintrack_hp;lintrack1];
           end
           if (max(linfield_hp) >= 3) % peak firing rate max larger than 3 Hz
               rm = [rm;linfield_hp'];
               pm = [pm;pos_hp'];
               tm = [tm;lintrack_hp'];
               cellidxm_hp = [cellidxm_hp; cind];
           end
      end
end
% PFC cell loop
for i = 1:ctxnum
      cind = ctxidx(i,:);
      if (length(linfields{day}{eprun})>= cind(1))
            if  (length(linfields{day}{eprun}{cind(1)})>= cind(2))
                linfield1 = linfields{day}{eprun}{cind(1)}{cind(2)};
            else 
                linfield1 =[];
            end
      else
            linfield1=[];
      end
      
      if ~isempty(linfield1)
           linfield_ctx = [];
           lintrack_ctx = [];
           pos_ctx = [];
           for track = 1:4
                temp1 = linfield1{track};
                pos1 = temp1(:,1);
                lintrack1 = ones(size(pos1))*track;
                occnormrate1 = temp1(:,5);
                linfield_ctx = [linfield_ctx;occnormrate1];
                pos_ctx = [pos_ctx;pos1];
                lintrack_ctx = [lintrack_ctx;lintrack1];
           end
           rm = [rm;linfield_ctx'];
           pm = [pm;pos_ctx'];
           tm = [tm;lintrack_ctx'];
           cellidxm_ctx = [cellidxm_ctx; cind];
      end
end

rm = rm'; %[nPosBin x nHPCells]
pm = pm';
tm = tm';
% remove reward-well regions, if wellcutoff > 0
for i = 1:4
    pm_traj = pm(find(tm == i));
    maxpos = max(max(pm_traj));
    rm(find(tm == i & pm <= wellcutoff)) = 0;
    rm(find(tm == i & pm >= maxpos-wellcutoff)) = 0;
end

rm(find(isnan(rm))) = 0 ; % exclude NaN
cellidxm  = [cellidxm_hp;cellidxm_ctx];

hpnum = length(cellidxm_hp(:,1)); % update cell number
ctxnum = length(cellidxm_ctx(:,1));
%%
%-------------------create reactivation matrix templates----------------------%
for tr = unique(tm)'
    colid = find(tm(:,1) == tr);
    rm_tr = rm(colid,:);
    crm_tr = corr(rm_tr);
    crm_tr = crm_tr - eye(size(crm_tr)); % remove auto-correlation
    crm_tr(find(isnan(crm_tr))) = 0 ;
    crm{tr} = crm_tr;
end
%%
%--- load spikes and ripple time---%
spikes = loaddatastruct(dir, animalprefix, 'spikes', day); % get spikes
% get ripple time
load(sprintf('%s%srippletime0%d.mat',dir,animalprefix,day));
rip = ripple{day}{ep}; 
riptimes(:,1) = rip.starttime;
riptimes(:,2) = rip.endtime;    
rip_starttime = 1000*riptimes(:,1);  % in ms


dur = 1000*(riptimes(:,2) - riptimes(:,1));
keepidx = find(dur >= 5*tBinSz);%at least 5 bins, 50 ms for 10ms bins; exclude events < 50ms
rip_starttime = rip_starttime(keepidx);
riptimes = riptimes(keepidx,:);
%%
% loop
if ~isempty(riptimes)
    celldata = [];
    spikecounts = [];
    % cell loop, measure active cells during each event
    for cellcount = 1:hpnum+ctxnum
        index = [day,ep,cellidxm(cellcount,:)] ;
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
        else
            spiketimes = [];
        end
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes);
            spikebins = spikebins(validspikes);
            tmpcelldata = [spiketimes spikebins];
        end
        if ~isempty(spiketimes)
            tmpcelldata(:,3) = cellcount;
        else 
            tmpcelldata = [0 0 cellcount];
        end
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikecounts = [spikecounts; spikecount];
    end
    cellcounts = sum((spikecounts(1:hpnum,:) > 0));
    ctxcellcounts = sum((spikecounts(hpnum+1:end,:) > 0));
    eventindex = find(cellcounts >= cellcountthresh & ctxcellcounts >= cellcountthresh); % at least 5 CA1 and 5 PFC cells active
    
    % event loop
    revent = 0;% reset count
    for event = 1:length(eventindex)
        spikecount_event = spikecounts(:,eventindex(event)).*tBinSz/1000/(riptimes(eventindex(event),2)-riptimes(eventindex(event),1));
        event
        cellsi = celldata(find(celldata(:,2)==eventindex(event)),3); 
        [cellsi,ia] = unique(cellsi,'first');
        [~,sortorder] = sort(ia);
        event_cellSeq = cellsi(sortorder);
        tmpind = find(celldata(:,2) == eventindex(event));
        spiketimes = celldata(tmpind,1);
        cellindex = celldata(tmpind,3);
        %-----create the event matrix during SWRs (spkT{cells}.spiketimes) -----%
        for cell = event_cellSeq'
            validspikeidx = find(cellindex == cell);
            spkT{cell} = spiketimes(validspikeidx).*1000;
        end     
        
        startevent = riptimes(eventindex(event),1).*1000;
        endevent = riptimes(eventindex(event),2).*1000;
        timebins = startevent:tBinSz:endevent; % timebins are the binedges
        nTBin = length(timebins)-1;
        nCell = hpnum+ctxnum;
        spkPerBin = zeros(1,nTBin, nCell); % keep the inactive cells as 0s.
        spkPerBin_raw = zeros(1,nTBin, nCell); % keep the inactive cells as 0s.
       
        for nn  = 1:hpnum+ctxnum
            cellInd = nn; %current cell
            if length(spkT) >= cellInd
                if ~isempty(spkT{cellInd})
%                     spkPerBin(1,:,cellInd) = histcounts(spkT{cellInd}, timebins); %[1 x nTBin x nCell]
                    temp = histc(spkT{cellInd}, timebins); %[1 x nTBin x nCell]
                    temp1 = smoothvect(temp(1:end-1), g1);
                    spkPerBin(1,:,cellInd) = temp1;
                    spkPerBin_raw(1,:,cellInd) = temp(1:end-1);
                end
            end
        end
        nSpkPerTBin = squeeze(sum(spkPerBin_raw,3)); %[nTBin x 1] number of spikes in tBin  
        nonzerobins = find(nSpkPerTBin > 0);
        
        %-----create the event matrix during SWRs-----%
        cswr = corr(squeeze(spkPerBin(1,:,cellsi)));
        cswr = cswr - eye(size(cswr)); % no auto-correlation
        cswr(find(isnan(cswr))) = 0 ; % exclude NaN
        % get active cell index
        active_hpid = find(cellsi <= hpnum);
        active_hpcellnum = length(active_hpid);
        active_ctxid = find(cellsi > hpnum);
        active_ctxcellnum = length(active_ctxid);
        
        % trajectory loop
        for tr = unique(tm)'
            temp = crm{tr};
            crm_ev = temp(cellsi,cellsi);
            crm_vec = [];
            cswr_vec = [];
            for i = 1:active_hpcellnum
                cswr_vec = [cswr_vec,cswr(i,1+active_hpcellnum:end)]; % SWR; only cross-regional correlations are used; within regional correlation not used
                crm_vec = [crm_vec, crm_ev(i,1+active_hpcellnum:end)];% RUN
            end
            R_rm_swr(tr) = corr(cswr_vec',crm_vec'); % r-value; R = C(RUN)* C(SWR); R is calculated for each trajectory
        end
        
        %-------Shuffling to get the pvalue for each traj------%
        scorr = [];
        % trajectory loop
        for tr = unique(tm)'
            temp = crm{tr};
            crm_ev = temp(cellsi,cellsi);
            for iteration = 1:1500 % 1500 shuffles, hard coded
                for nn = cellsi'
                    permbins = randperm(length(nSpkPerTBin)); % randomly shuffle spike time during the SWR
                    temp = squeeze(spkPerBin_raw(1,:,nn));
                    temp_shuffle = temp(permbins);
                    temp1 = smoothvect(temp_shuffle, g1);
                    tmpspkPerBin(1,:,nn) = temp1;
                end
                % create C(SWR) based on shuffled data
                scswr = corr(squeeze(tmpspkPerBin(1,:,cellsi)));
                scswr = scswr - eye(size(scswr));
                scswr(find(isnan(scswr))) = 0 ; % exclude NaN if exist
                
                % calculate r-value based on shuffled data
                scswr_vec = [];
                crm_vec = [];
                for i = 1:active_hpcellnum
                    scswr_vec = [scswr_vec,scswr(i,1+active_hpcellnum:end)];
                    crm_vec = [crm_vec, crm_ev(i,1+active_hpcellnum:end)];
                end
                sR_rm_swr(tr,iteration) = corr(scswr_vec',crm_vec');% shuffled r-value
                clear tmpspkPerBin
            end
            
            % calculate p-value
            if R_rm_swr(tr) >= 0 
                pvalue(tr) =sum(R_rm_swr(tr) < sR_rm_swr(tr,:))/length(sR_rm_swr(tr,:)); 
            else
                pvalue(tr) =sum(R_rm_swr(tr) > sR_rm_swr(tr,:))/length(sR_rm_swr(tr,:)); 
            end
        end
        [minP,tidx] = min(pvalue);
        % decoded trajectory is the one with minimum pvalue
        if minP < 0.05
            decode_traj = tidx; % significant trajectory
        else
            decode_traj = 0;% no significant trajectory
        end
        
        % structure result
        reactivationtraj{day}{ep}.eventinfo{event}.eventime = [startevent, endevent]./1000;
        reactivationtraj{day}{ep}.eventinfo{event}.activecell_id = cellsi;
        reactivationtraj{day}{ep}.eventinfo{event}.activecell_hpid = cellsi(active_hpid);
        reactivationtraj{day}{ep}.eventinfo{event}.activecell_ctxid = cellsi(active_ctxid);
        reactivationtraj{day}{ep}.eventinfo{event}.totalcell = nCell;
        reactivationtraj{day}{ep}.eventinfo{event}.activecell_info = cellidxm(cellsi,:); 
        reactivationtraj{day}{ep}.eventinfo{event}.CorrMat_SWR = cswr;
        reactivationtraj{day}{ep}.eventinfo{event}.CorrMat_Behav = crm;
        reactivationtraj{day}{ep}.eventinfo{event}.decodedtraj = decode_traj;
        reactivationtraj{day}{ep}.eventinfo{event}.rvalue_all = R_rm_swr;
        reactivationtraj{day}{ep}.eventinfo{event}.pvalue_shuffle_all = pvalue;
        reactivationtraj{day}{ep}.eventinfo{event}.reactivationevent = sign(decode_traj);
        if reactivationtraj{day}{ep}.eventinfo{event}.reactivationevent % significant event
            reactivationtraj{day}{ep}.eventinfo{event}.rvalue = R_rm_swr(decode_traj);
            revent = revent + 1; % increase count
        end
    end
    reactivationtraj{day}{ep}.candeventnum = event;
    reactivationtraj{day}{ep}.sigeventnum = revent;
    reactivationtraj{day}{ep}.sigeventprc = revent./event;
end

%%
%---save date ---%
if savedata
   save(sprintf('%s%sreactivationtraj_hpctx_%02d.mat', savedir,animalprefix,day), 'reactivationtraj');
end 
        
        
        


