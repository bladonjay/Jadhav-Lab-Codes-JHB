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
    trialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10];
    trialmat(trialmat(:,4)==0,:)=[]; %remove incorrect trials
    [temp,sortinds]=sortrows(trialmat,[3 1]); cutspot=find(diff(trialmat(sortinds,3)~=0)); odorid=trialmat(:,3);
    % okay now get the perievent raster, lets do the buzsaki way thats a
    % 100 msec smoothing kernel for each trial and grab -2 seconds to +2
    % seconds and patch the mean curve
    
    for i=1:length(SuperRat(ses).units)
        
        % get spikes
        [~,~,~,~,~,spikets]=event_spikes(SuperRat(ses).units(i).ts,trialmat(:,1),abs(timeedges(1)),timeedges(end));
        % now smooth each into a vector
        curves=cellfun(@(a) SmoothMat2(histcounts(a,timeedges),[1000 0],100), spikets, 'UniformOutput', false);
        
        % or unsmoothed
        %curves=cellfun(@(a)histcounts(a,-2:.001:2), spikets, 'UniformOutput', false);
        curvemat=cell2mat(curves');
        
        % now get the trial rates for object selectivity
        [~,trspikes]=event_spikes(SuperRat(ses).units(i).ts,trialmat(:,1),0,trialmat(:,2)-trialmat(:,1));
        % and now get the pre trial rates for a trial-length matched control
        [~,pretrspikes]=event_spikes(SuperRat(ses).units(i).ts,trialmat(:,1)-(trialmat(:,2)-trialmat(:,1)),0,trialmat(:,2)-trialmat(:,1));
        
        coderp=ranksum(trspikes(trialmat(:,3)==1),trspikes(trialmat(:,3)==0));
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
                mean(trialmat(:,2)-trialmat(:,1)),responderp,coderp));
        end
    end
end

%%
% this is claires way
useblocks=false;
bootct=500;

% this adds the following fields
% OdorRates: a nx3 matrix, 1 rate, 2 nspikes, 3 lr10 (doesnt get you the
% CI10
% odorSelective: a 4 element vector, 1 diff in mean rates, 2 pval, 3 is
% p<05
% odorMeans: mean rate for each odor
% odorResponsive: rate at odor, mean chg from before, std chg from before,
% pval


for ses=1:length(SuperRat)
    tic
    % so maybe something like a left plot mean rate patch, then right plot
    % a colormap of abs diff/sum?
    trialdata=SuperRat(ses).trialdata;
    % first grab the tuning curves during the one second delay
    % mat is 1 start, 2 end 3 left right 4 correct incorrect and 5 epoch
    trialmat=[trialdata.sniffstart trialdata.sniffend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
    goodtrials=trialmat(:,4)==1; odorid=trialmat(:,3);
    % okay now get the perievent raster, lets do the buzsaki way thats a
    % 100 msec smoothing kernel for each trial and grab -2 seconds to +2
    % seconds and patch the mean curve
    wb=waitbar(0,'Starting to run cells');
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorRates'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorSelective'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorMeans'); end
    try, SuperRat(ses).units=rmfield(SuperRat(ses).units,'OdorResponsive'); end
    
    for i=1:length(SuperRat(ses).units)
        %build  the spikes
        [totalspikes,spkevs,~,spikeinds]=event_spikes(SuperRat(ses).units(i).ts(:,1),trialmat(:,1),0,trialmat(:,2)-trialmat(:,1));
        [~,prespkevs]=event_spikes(SuperRat(ses).units(i).ts(:,1),trialmat(:,1)-(trialmat(:,2)-trialmat(:,1)),0,trialmat(:,2)-trialmat(:,1));
        spkcts=cellfun(@(a) length(a), spikeinds);
        %curves=cellfun(@(a) SmoothMat2(histcounts(a,0:.001:1),[1000 0],100), controls, 'UniformOutput', false);
        %curvemat=cell2mat(curves');
        if useblocks
            epochcounts=unique(trialmat(:,5)); % how many sessions are in this recording day?
            SIdata=[]; rlmeans=[];
            for ep=1:length(unique(epochcounts)) % calc for each session
                % get all the spikes for the correct trials and this epoch
                goodtrials=trialmat(:,4)==1 & trialmat(:,5)==epochcounts(ep);
                if totalspikes>=length(spkevs) % if the cell has more spikes than there are trials
                    % get the mean rates for each
                    % LR10 left first, right second
                    rlmeans(ep,1)=nanmean(spkevs(goodtrials & trialmat(:,3)==1)); % and this is spikes per second
                    rlmeans(ep,2)=nanmean(spkevs(goodtrials & trialmat(:,3)==0));
                    % and an effect size just to see
                    %rlmeans(ep,3)=dprime(spkevs(usetrials & trialmat(:,3)==1),spkevs(usetrials & trialmat(:,3)==0));
                    % and selectivity index
                    [SIdata(ep,1),SIdata(ep,2)]=SelectivityIndex(spkevs(goodtrials),trialmat(goodtrials,3),bootct); % 500 boots
                    SIdata(ep,3)=SIdata(ep,2) <= 0.05 & ~isnan(SIdata(ep,1)); % how many boots
                    
                    % and the responsivity index (and i'll do a 'change
                    % from baseline'
                    
                else
                    rlmeans(ep,1:3)=nan; SIdata(ep,1:3)=nan; % if not enouch spikes, its nanned out
                end
                
                % is the cell odor responsive
                pokeresponse=nanmean(spksevs(goodtrials));
                pokeresponse(2)=nanmean(spkevs(goodtrials)-prespkevs(goodtrials)); % + if elevated, - if depressed
                pokeresponse(3)=nanstd(spkevs(goodtrials)-prespkevs(goodtrials)); % want to know if its consistent
                pokeresponse=signrank(spkevs(goodtrials),prespkevs(goodtrials));
                SuperRat(ses).units(i).OdorResponsive=[pokeresponse pokeresponse];
                
                SuperRat(ses).units(i).OdorRates{ep}=[spkevs(goodtrials) spkcts(goodtrials)' trialmat(goodtrials,3)];
                SuperRat(ses).units(i).OdorMeans(ep,:)=rlmeans; % upload to the struct
                SuperRat(ses).units(i).OdorSelective(ep,:)=SIdata;
            end
        else
            if totalspikes>=length(spkevs) % if the cell spikes on at least half the trials
                % get the mean rates for each
                % LR10 left first, right second
                rlmeans=nanmean(spkevs(goodtrials & odorid==1)); % and this is spikes per second
                rlmeans(1,2)=nanmean(spkevs(goodtrials & odorid==0));
                % and an effect size just to see
                SIdata=abs(diff(fliplr(rlmeans)))/sum(rlmeans); % flip cause i want - if 2 is larger
                % and selectivity index
                [~,SIdata(1,2)]=SelectivityIndex(spkevs(goodtrials),trialmat(goodtrials,3),bootct); % 500 boots
                SIdata(1,3)=SIdata(2) <= 0.05 & ~isnan(SIdata(1)); % how many boots
                SIdata(1,4)=dprime(spkevs(goodtrials & trialmat(:,3)==1),spkevs(goodtrials & trialmat(:,3)==0));
            else
                rlmeans(1,1:3)=nan; SIdata(1,1:4)=nan; % if not enouch spikes, its nanned out
            end
            
            % poke responsiveness= should it be paired or not?
            % should probably do this for each of the odors
            pokeresponse=nanmean(spkevs(goodtrials & odorid==1)); % what is the sampling rate
            pokeresponse(1,2)=nanmean(spkevs(goodtrials & odorid==1)-prespkevs(goodtrials & odorid==1)); % + if elevated, - if depressed
            pokeresponse(1,3)=dprime(spkevs(goodtrials & odorid==1),prespkevs(goodtrials & odorid==1)); % want to know if its consistent (normalized effect size)
            pokeresponse(1,4)=signrank(spkevs(goodtrials & odorid==1),prespkevs(goodtrials & odorid==1)); % is it significant
            pokeresponse(2,1)=nanmean(spkevs(goodtrials & odorid==0)); % what is the sampling rate
            pokeresponse(2,2)=nanmean(spkevs(goodtrials & odorid==0)-prespkevs(goodtrials & odorid==0)); % + if elevated, - if depressed
            pokeresponse(2,3)=dprime(spkevs(goodtrials & odorid==0),prespkevs(goodtrials & odorid==0)); % want to know if its consistent
            pokeresponse(2,4)=signrank(spkevs(goodtrials & odorid==0),prespkevs(goodtrials & odorid==0)); % is it significant
            SuperRat(ses).units(i).OdorResponsive=pokeresponse;
            
            
            SuperRat(ses).units(i).OdorRates=[spkevs(goodtrials)' spkcts(goodtrials)' trialmat(goodtrials,3)];
            SuperRat(ses).units(i).OdorMeans=rlmeans; % upload to the struct
            SuperRat(ses).units(i).OdorSelective=SIdata;
        end
        
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
        trialmat=[trialdata.rewardstart trialdata.rewardend trialdata.leftright10 trialdata.CorrIncorr10 trialdata.EpochInds(:,2)];
        [temp,sortinds]=sortrows(trialmat,[3 1]); odorid=trialmat(:,3);
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
            [~,spkevs,~,~,~,controls]=event_spikes(SuperRat(ses).units(i).ts(:,1),trialmat(:,1),0,trialmat(:,2)-trialmat(:,1));
            %curves=cellfun(@(a) SmoothMat2(histcounts(a,0:.001:1),[1000 0],100), controls, 'UniformOutput', false);
            %curvemat=cell2mat(curves');
            if useblocks
                epochcounts=unique(trialmat(:,5)); % how many sessions are in this recording day?
                SIdata=[]; rlmeans=[];
                for ep=1:length(unique(epochcounts)) % calc for each session
                    % get all the spikes for the correct trials and this epoch
                    goodtrials=trialmat(:,4)==1 & trialmat(:,5)==epochcounts(ep); %
                    if nanmean(any(spkevs(goodtrials)))>=0.5 % if the cell spikes on at least half the trials
                        % get the mean rates for each
                        % LR10 left first, right second
                        rlmeans(ep,1)=nanmean(spkevs(goodtrials & trialmat(:,3)==1)); % and this is spikes per second
                        rlmeans(ep,2)=nanmean(spkevs(goodtrials & trialmat(:,3)==0));
                        % and an effect size just to see
                        %rlmeans(ep,3)=dprime(spkevs(usetrials & trialmat(:,3)==1),spkevs(usetrials & trialmat(:,3)==0));
                        % and selectivity index
                        [SIdata(ep,1),SIdata(ep,2)]=SelectivityIndex(spkevs(goodtrials),trialmat(goodtrials,3),bootct); % 500 boots
                        SIdata(ep,3)=SIdata(ep,2) <= 0.05 & ~isnan(SIdata(ep,1)); % how many boots
                    else
                        rlmeans(ep,1:3)=nan; SIdata(ep,1:3)=nan; % if not enouch spikes, its nanned out
                    end
                    SuperRat(ses).units(i).RewardRates{ep}=[spkevs(goodtrials) trialmat(goodtrials,3)];
                    SuperRat(ses).units(i).RewardMeans(ep,:)=rlmeans; % upload to the struct
                    SuperRat(ses).units(i).RewardSelective(ep,:)=SIdata;
                end
            else
                goodtrials=trialmat(:,4)==1;  % only correct trials
                if nanmean(any(spkevs))>=0.5 % if the cell spikes on at least half the trials
                    % get the mean rates for each
                    % LR10 left first, right second
                    rlmeans=nanmean(spkevs(goodtrials & trialmat(:,3)==1)); % and this is spikes per second
                    rlmeans(1,2)=nanmean(spkevs(goodtrials & trialmat(:,3)==0));
                    % and an effect size just to see
                    %rlmeans(ep,3)=dprime(spkevs(usetrials & trialmat(:,3)==1),spkevs(usetrials & trialmat(:,3)==0));
                    % and selectivity index
                    [SIdata,SIdata(1,2)]=SelectivityIndex(spkevs(goodtrials),trialmat(goodtrials,3),bootct); % 500 boots
                    SIdata(1,3)=SIdata(2) <= 0.05 & ~isnan(SIdata(1)); % how many boots
                else
                    rlmeans(1,1:3)=nan; SIdata(1,1:3)=nan; % if not enouch spikes, its nanned out
                end
                SuperRat(ses).units(i).RewardRates=[spkevs(goodtrials)' trialmat((goodtrials),3)];
                SuperRat(ses).units(i).RewardMeans=rlmeans; % upload to the struct
                SuperRat(ses).units(i).RewardSelective=SIdata;
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
