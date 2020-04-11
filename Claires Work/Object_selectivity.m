% now lets add some object selectivity onto this

% first, how many cells show object selectivity, lets see if we have the
% data in the behavior first


%load('E:\Claire Data\ClaireDataFull-LinPos-LongTrack.mat')


% object epochs are 'sniff' and the right and left 'correct' are the sides,
% but if its inco0rrect, you have to flip the odor
%% will ahve to use cs33 & cs 34 because 44 doesn thave odor identities\



%for ses=1:length(SuperRat)
for ses=1:length(SuperRat)
    
    timeedges=-2:.001:2; timebins=-2:.001:1.999;
    % so maybe something like a left plot mean rate patch, then right plot
    % a colormap of abs diff/sum?
    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    mytrialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10];
    mytrialmat(mytrialmat(:,4)==0,:)=[]; %remove incorrect trials
    [temp,sortinds]=sortrows(mytrialmat,[3 1]); cutspot=find(diff(mytrialmat(sortinds,3)~=0)); odorid=mytrialmat(:,3);
    % okay now get the perievent raster, lets do the buzsaki way thats a
    % 100 msec smoothing kernel for each trial and grab -2 seconds to +2
    % seconds and patch the mean curve
    
    for i=1:length(SuperRat(ses).units)
        
        % get spikes
        [~,~,~,~,~,spikets]=event_spikes(SuperRat(ses).units(i).ts,mytrialmat(:,1),abs(timeedges(1)),timeedges(end));
        % now smooth each into a vector
        curves=cellfun(@(a) SmoothMat2(histcounts(a,timeedges),[1000 0],100), spikets, 'UniformOutput', false);
        
        % or unsmoothed
        %curves=cellfun(@(a)histcounts(a,-2:.001:2), spikets, 'UniformOutput', false);
        curvemat=cell2mat(curves');
        
        % now get the trial rates for object selectivity
        [~,trspikes]=event_spikes(SuperRat(ses).units(i).ts,mytrialmat(:,1),0,mytrialmat(:,2)-mytrialmat(:,1));
        % and now get the pre trial rates for a trial-length matched control
        [~,pretrspikes]=event_spikes(SuperRat(ses).units(i).ts,mytrialmat(:,1)-(mytrialmat(:,2)-mytrialmat(:,1)),0,mytrialmat(:,2)-mytrialmat(:,1));
        
        coderp=ranksum(trspikes(mytrialmat(:,3)==1),trspikes(mytrialmat(:,3)==0));
        responderp=signrank(trspikes,pretrspikes);
        if coderp<.05 || responderp<.05
            figure;
            sp(1)=subplot(2,1,1);
            imagesc(timebins, 1:size(curves,2), curvemat(sortinds,:)); hold on;
            plot([timebins(1) timebins(end)],[cutspot cutspot],'k','LineWidth',2);
            sp(2)=subplot(2,1,2);
            plot(timebins,nanmean(curvemat(odorid==1,:)),'r','LineWidth',3);
            
            patch([timebins fliplr(timebins)],[nanmean(curvemat(odorid==1,:))+nanstd(curvemat(odorid==1,:),0,1)...
                fliplr(nanmean(curvemat(odorid==1,:))-nanstd(curvemat(odorid==1,:),0,1))],[.85 .325 .098],...
                'LineStyle','none','FaceAlpha',.5)
            hold on; plot(timebins,nanmean(curvemat(odorid==0,:)),'b','LineWidth',3);
            patch([timebins fliplr(timebins)],[nanmean(curvemat(odorid==0,:))+nanstd(curvemat(odorid==0,:),0,1)...
                fliplr(nanmean(curvemat(odorid==0,:))-nanstd(curvemat(odorid==0,:),0,1))],[0    0.4470    0.7410],...
                'LineStyle','none','FaceAlpha',.5)
            linkaxes(sp,'x');
            % now get the overall rate
            title(sprintf('Mean Sampling period = %.2f P for responsive is %.2f p for selective is %.2f',...
                mean(mytrialmat(:,2)-mytrialmat(:,1)),responderp,coderp));
        end
    end
end

%%
% this is claires way
useblocks=false;
bootct=1000;

% this adds the following fields:


% OdorRates: a nx3 matrix, 1 rate, 2 nspikes, 3 lr10 (doesnt get you the CI10
% odorSelective: a 4 element vector, 1 diff in mean rates, 2 pval, 3 is
% p<05, 4 is the dprime effect size
% odorMeans: mean rate for each odor
% odorResponsive: rate at odor, mean chg from before, std chg from before, pval


for ses=1:length(SuperRat)
    tic

    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    % mat is 1 start, 2 end 3 left right 4 correct incorrect and 5 epoch
    fulltrialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
    trialcorrect=fulltrialmat(:,4)==1 & ismember(fulltrialmat(:,5),SuperRat(ses).RunEpochs);
    trialmat=fulltrialmat(trialcorrect,:); % only take correct trials
    

    wb=waitbar(0,'Starting to run cells');
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorRates'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorSelective'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorMeans'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorResponsive'); end
    
    for i=1:length(SuperRat(ses).units)
        % first grab all the events and spikes
        [~,spkevs,~,trspks]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
            trialmat(:,1),0,trialmat(:,2)-trialmat(:,1));
        % remove trials from blocks where the cell has fewer spikes than
        % trials
        spknums=cellfun(@(a) length(a), trspks);
        spikesperblock=accumarray(trialmat(:,5),spknums);
        trialsperblock=accumarray(trialmat(:,5),1);
        %find those inds and remove them
        keepblocktrials=sum(find(spikesperblock>trialsperblock/2)'==trialmat(:,5),2)>0;
        mytrialmat=trialmat(keepblocktrials,:);
        myspkevs=spkevs(keepblocktrials);
        % preallocate
        pokeresponse=nan(2,4); rlmeans=zeros(1,2); Selectivitydata=nan(1,4);
        if ~isempty(mytrialmat) % cell has to have more spikes than there were trials
            % now split out by odor
            odorid=mytrialmat(:,3);
            % get the spike rate vectors
            [totalspikes,spkevs,~,spikeinds]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
                mytrialmat(:,1),0,mytrialmat(:,2)-mytrialmat(:,1));
            [~,prespkevs]=event_spikes(SuperRat(ses).units(i).ts(:,1),...
                mytrialmat(:,1)-(mytrialmat(:,2)-mytrialmat(:,1)),0,mytrialmat(:,2)-mytrialmat(:,1));
            spkcts=cellfun(@(a) length(a), spikeinds);

            % get the mean rates for each
            % LR10 left first, right second
            rlmeans=nanmean(spkevs(odorid==1)); % and this is spikes per second
            rlmeans(1,2)=nanmean(spkevs(odorid==0));
            % and an effect size just to see
            Selectivitydata=abs(diff(fliplr(rlmeans)))/sum(rlmeans); % flip cause i want - if 2 is larger
            % and selectivity index
            [~,Selectivitydata(1,2)]=SelectivityIndex(spkevs,mytrialmat(:,3),bootct); % 500 boots
            Selectivitydata(1,3)=Selectivitydata(2) <= 0.05 & ~isnan(Selectivitydata(1)); % how many boots
            Selectivitydata(1,4)=dprime(spkevs(mytrialmat(:,3)==1),spkevs(mytrialmat(:,3)==0));
            
            % sign rank for each odor separately, then take the best
            % answer
            for r=1:2
                pokeresponse((3-r),1)=nanmean(spkevs(odorid==(2-r))); % what is the sampling rate
                pokeresponse((3-r),2)=nanmean(spkevs(odorid==(2-r))-prespkevs(odorid==(2-r))); % + if elevated, - if depressed
                pokeresponse((3-r),3)=dprime(spkevs(odorid==(2-r)),prespkevs(odorid==(2-r))); % want to know if its consistent (normalized effect size)
                pokeresponse((3-r),4)=signrank(spkevs(odorid==(2-r)),prespkevs(odorid==(2-r))); % is it significant
            end
        end
        SuperRat(ses).units(i).OdorResponsive=pokeresponse;
        SuperRat(ses).units(i).OdorRates=[myspkevs' mytrialmat(:,3)];
        SuperRat(ses).units(i).OdorMeans=rlmeans; % upload to the struct
        SuperRat(ses).units(i).OdorSelective=Selectivitydata;
        
        waitbar(i/length(SuperRat(ses).units),wb,sprintf('Running unit %d',i));
    end
    close(wb);
    % record how long it took (it usually is super short)
    fprintf('Ses %d took %.2f seconds \n',ses,toc);
    
end

%%now plot the spatial characteristics of those cells


%% Repeat for reward poke firing

%************ DONT HAVE ALL THE REWARD DATA THOUGH ****************

% this is claires way
useblocks=false;
bootct=500;

for ses=1:length(SuperRat)
    tic
    % so maybe something like a left plot mean rate patch, then right plot
    % a colormap of abs diff/sum?
    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    % mat is 1 start, 2 end 3 left right 4 correct incorrect and 5 epoch\
    if isfield(trialdata,'rewardstart')
        mytrialmat=[trialdata.rewardstart trialdata.rewardend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
        [temp,sortinds]=sortrows(mytrialmat,[3 1]); odorid=mytrialmat(:,3);
        % okay now get the perievent raster, lets do the buzsaki way thats a
        % 100 msec smoothing kernel for each trial and grab -2 seconds to +2
        % seconds and patch the mean curve
        rlmeans=nan(length(SuperRat(ses).units),2);
        wb=waitbar(0,'Starting to run cells');
        try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'RewardRates'); end
        try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'RewardMeans'); end
        try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'RewardSelective'); end
        
        for i=1:length(SuperRat(ses).units)
            
            %build  the spikes
            [~,spkevs,~,~,~,controls]=event_spikes(SuperRat(ses).units(i).ts(:,1),mytrialmat(:,1),0,mytrialmat(:,2)-mytrialmat(:,1));
            %curves=cellfun(@(a) SmoothMat2(histcounts(a,0:.001:1),[1000 0],100), controls, 'UniformOutput', false);
            %curvemat=cell2mat(curves');
            if useblocks
                epochcounts=unique(mytrialmat(:,5)); % how many sessions are in this recording day?
                Selectivitydata=[]; rlmeans=[];
                for ep=1:length(unique(epochcounts)) % calc for each session
                    % get all the spikes for the correct trials and this epoch
                    goodtrials=mytrialmat(:,4)==1 & mytrialmat(:,5)==epochcounts(ep); %
                    if nanmean(any(spkevs(goodtrials)))>=0.5 % if the cell spikes on at least half the trials
                        % get the mean rates for each
                        % LR10 left first, right second
                        rlmeans(ep,1)=nanmean(spkevs(goodtrials & mytrialmat(:,3)==1)); % and this is spikes per second
                        rlmeans(ep,2)=nanmean(spkevs(goodtrials & mytrialmat(:,3)==0));
                        % and an effect size just to see
                        %rlmeans(ep,3)=dprime(spkevs(usetrials & trialmat(:,3)==1),spkevs(usetrials & trialmat(:,3)==0));
                        % and selectivity index
                        [Selectivitydata(ep,1),Selectivitydata(ep,2)]=SelectivityIndex(spkevs(goodtrials),mytrialmat(goodtrials,3),bootct); % 500 boots
                        Selectivitydata(ep,3)=Selectivitydata(ep,2) <= 0.05 & ~isnan(Selectivitydata(ep,1)); % how many boots
                    else
                        rlmeans(ep,1:3)=nan; Selectivitydata(ep,1:3)=nan; % if not enouch spikes, its nanned out
                    end
                    SuperRat(ses).units(i).RewardRates{ep}=[spkevs(goodtrials) mytrialmat(goodtrials,3)];
                    SuperRat(ses).units(i).RewardMeans(ep,:)=rlmeans; % upload to the struct
                    SuperRat(ses).units(i).RewardSelective(ep,:)=Selectivitydata;
                end
            else
                goodtrials=mytrialmat(:,4)==1;  % only correct trials
                if nanmean(any(spkevs))>=0.5 % if the cell spikes on at least half the trials
                    % get the mean rates for each
                    % LR10 left first, right second
                    rlmeans=nanmean(spkevs(goodtrials & mytrialmat(:,3)==1)); % and this is spikes per second
                    rlmeans(1,2)=nanmean(spkevs(goodtrials & mytrialmat(:,3)==0));
                    % and an effect size just to see
                    %rlmeans(ep,3)=dprime(spkevs(usetrials & trialmat(:,3)==1),spkevs(usetrials & trialmat(:,3)==0));
                    % and selectivity index
                    [Selectivitydata,Selectivitydata(1,2)]=SelectivityIndex(spkevs(goodtrials),mytrialmat(goodtrials,3),bootct); % 500 boots
                    Selectivitydata(1,3)=Selectivitydata(2) <= 0.05 & ~isnan(Selectivitydata(1)); % how many boots
                else
                    rlmeans(1,1:3)=nan; Selectivitydata(1,1:3)=nan; % if not enouch spikes, its nanned out
                end
                SuperRat(ses).units(i).RewardRates=[spkevs(goodtrials)' mytrialmat((goodtrials),3)];
                SuperRat(ses).units(i).RewardMeans=rlmeans; % upload to the struct
                SuperRat(ses).units(i).RewardSelective=Selectivitydata;
            end
            
            waitbar(i/length(SuperRat(ses).units),wb,sprintf('Running unit %d',i));
        end
        close(wb);
        % record how long it took (it usually is super short)
        fprintf('Ses %d took %.2f seconds \n',ses,toc);
    else
        fprintf('Ses %d doesnt have reward events \n',ses);
        SuperRat(ses).units(:).RewardRates=nan(1,3);
        
    end
    
end

%% lets get a list of each cells session and name

SuperUnits=orderfields(SuperRat(1).units);
%SuperUnits=rmfield(SuperUnits,'csi');
thisdate=repmat({[SuperRat(1).name ' day ' num2str(SuperRat(1).daynum)]},length(SuperUnits),1);

[SuperUnits.sess]=thisdate{:};
for i=2:length(SuperRat)
    theseunits=orderfields(SuperRat(i).units);
    % for some reason csi is in some of these structs
    if isfield(theseunits,'csi'), theseunits=rmfield(theseunits,'csi'); end
    thisdate=repmat({[SuperRat(i).name ' day ' num2str(SuperRat(i).daynum)]},length(theseunits),1);
    [theseunits.sess]=thisdate{:};
    SuperUnits=[SuperUnits theseunits];
end

% and the comparison to claires dataset
clairedata=load('E:\Claire Data\ClaireObjResults.mat');
clairedata.names={'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
% so i have to remove all sessions with rat numbers
%  6 and 7

% okay now need to track down all the cell she thinks code...
clear ClaireRespCA1;

for i=1:length(clairedata.npCellsCA1)
    ClaireRespCA1{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.npCellsCA1(i,1)},...
        clairedata.npCellsCA1(i,2),clairedata.npCellsCA1(i,3),...
        clairedata.npCellsCA1(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.npCellsCA1(i,1)},clairedata.npCellsCA1(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.npCellsCA1(i,3) & [SuperUnits.unitnum]==clairedata.npCellsCA1(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireRespCA1{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).OdorResponsive;
        [~,useind]=min(odordata(:,4));
        ClaireRespCA1{i,3}=odordata(useind,:);
        ClaireRespCA1{i,4}=odordata(useind,4)<.05;
        
    elseif clairedata.npCellsCA1(i,1)==6 || clairedata.npCellsCA1(i,1)==7
        ClaireRespCA1{i,2}='didnt analyze session';
    else
        ClaireRespCA1{i,2}='couldnt find cell';
    end
end


%okay now need to track down all the cell she thinks code...
clear ClaireSelCA1;
for i=1:length(clairedata.selectiveCA1)
    ClaireSelCA1{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.selectiveCA1(i,1)},...
        clairedata.selectiveCA1(i,2),clairedata.selectiveCA1(i,3),...
        clairedata.selectiveCA1(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.selectiveCA1(i,1)},clairedata.selectiveCA1(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.selectiveCA1(i,3) & [SuperUnits.unitnum]==clairedata.selectiveCA1(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireSelCA1{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).OdorResponsive;
        [~,useind]=min(odordata(:,4));
        ClaireSelCA1{i,3}=odordata(useind,:);
        ClaireSelCA1{i,4}=odordata(useind,4)<.05;
        
    elseif clairedata.selectiveCA1(i,1)==6 || clairedata.selectiveCA1(i,1)==7
        ClaireSelCA1{i,2}='didnt analyze session';
    else
        ClaireSelCA1{i,2}='couldnt find cell';
    end
end


% okay now need to track down all the cell she thinks code...
clear ClaireRespPFC;
for i=1:length(clairedata.npCellsPFC)
    ClaireRespPFC{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.npCellsPFC(i,1)},...
        clairedata.npCellsPFC(i,2),clairedata.npCellsPFC(i,3),...
        clairedata.npCellsPFC(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.npCellsPFC(i,1)},clairedata.npCellsPFC(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.npCellsPFC(i,3) & [SuperUnits.unitnum]==clairedata.npCellsPFC(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireRespPFC{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).OdorResponsive;
        [~,useind]=min(odordata(:,4));
        ClaireRespPFC{i,3}=odordata(useind,:);
        ClaireRespPFC{i,4}=odordata(useind,4)<.05;
        
    elseif clairedata.npCellsPFC(i,1)==6 || clairedata.npCellsPFC(i,1)==7
        ClaireRespPFC{i,2}='didnt analyze session';
    else
        ClaireRespPFC{i,2}='couldnt find cell';
    end
end


% okay now need to track down all the cell she thinks code...
clear ClaireSelPFC;
for i=1:length(clairedata.selectivePFC)
    ClaireSelPFC{i,1}=sprintf('%s ses %d tet %d cluster %d',...
        clairedata.names{clairedata.selectivePFC(i,1)},...
        clairedata.selectivePFC(i,2),clairedata.selectivePFC(i,3),...
        clairedata.selectivePFC(i,4));
    matchsess=contains({SuperUnits.sess},sprintf('%s day %d',...
        clairedata.names{clairedata.selectivePFC(i,1)},clairedata.selectivePFC(i,2)));
    matchcell=[SuperUnits.tet]==clairedata.selectivePFC(i,3) & [SuperUnits.unitnum]==clairedata.selectivePFC(i,4);
    matchind=find(matchsess & matchcell);
    if ~isempty(matchind)
        ClaireSelPFC{i,2}=sprintf('%s tet %d unit %d',SuperUnits(matchind).sess,...
            SuperUnits(matchind).tet, SuperUnits(matchind).unitnum);
        odordata=SuperUnits(matchind).OdorResponsive;
        [~,useind]=min(odordata(:,4));
        ClaireSelPFC{i,3}=odordata(useind,:);
        ClaireSelPFC{i,4}=odordata(useind,4)<.05;
        
    elseif clairedata.selectivePFC(i,1)==6 || clairedata.selectivePFC(i,1)==7
        ClaireSelPFC{i,2}='didnt analyze session';
    else
        ClaireSelPFC{i,2}='couldnt find cell';
    end
end





