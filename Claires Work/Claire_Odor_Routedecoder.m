%% so now the decoder analysis
%{
thhe gist is that i want to encode the odor sampling period, and then use
it to decode the trajectory identity at each step along each run.

% i can think of two ways to do this. first, run the exact code wenbo used
to decode ripples but in reverse.  the other way is to train on the odor
period and match 800 msec bins for each and every 'run' epoch, find the
most well predicted bin in that epoch, classify, and mark its position...
or mark every 800 msec bin, and aggregate the mean positions for each
decode.

% the rub i anticipate for the wenbo method is i'm not sure how to
calculate p(spike|odor) if the rates excede 1.


%}
%CellMatrix is a trial by cell matrix of nspikes
% trialIds is a column of left and right odors


% so this will be done for each PFC and dHPC, will ahve to match for cell
% counts, and will have to probably aggregate per area.  If need be i can
% do a rigotti style analyisis where i randomly pull from any animal, but
% i'd rather not do that




region={'CA1','PFC'};


% set random stream so results are always the same
rstream = RandStream('mt19937ar','Seed',16);
RandStream.setGlobalStream(rstream);
%% this is claires algorithm...

% we are not using this for the decoder right now

%{

warning('off','all');
winsize=1;
verbose=1;
nfolds = 5;
decoderprobs=[];
for ses=1:length(SuperRat)
    
    % we're going to predict the action of the animal, so correct and
    % incorrect trials are used
    odorOK= ismember(SuperRat(ses).trialdata.EpochInds(:,2),SuperRat(ses).RunEpochs);
    odorflip=~logical(SuperRat(ses).trialdata.CorrIncorr10(odorOK)); % wrong trials
    odorStarts=SuperRat(ses).trialdata.sniffstart(odorOK); % start time
    odorIDs=SuperRat(ses).trialdata.leftright10(odorOK); % odor id
    badtrials=isnan(odorIDs);
    % somehow nans snuck into the odor ids? wtf?
    odorStarts(badtrials)=[]; odorIDs(badtrials)=[]; odorflip(badtrials)=[];
    odorIDs(odorflip)=double(~odorIDs(odorflip)); % flip the sign of the wrong trials
    %
    
    odorSpkMat=[]; % grab big old matrix
    for j=1:length(SuperRat(ses).units)
        [~,odorSpkMat(:,j)]=event_spikes(SuperRat(ses).units(j).ts,...
            odorStarts,0,winsize);
    end
    badunits=sum(odorSpkMat>0,1)<10 | nanmean(odorSpkMat)>9; % less than ten active trials is bad
    % now index columns that are members of PFC or dHPC
    
    for i=1:length(region)
        inRegion=cellfun(@(a) contains(a,region{i},'IgnoreCase',true), {SuperRat(ses).units.area});
        SpkMat=odorSpkMat(:,inRegion & ~badunits).*winsize; % p(spike| 1 ms timebin)
        
        % true decoding
        cv = cvpartition(length(odorIDs), 'kfold',nfolds);
        clear fract_correct;
        for k=1:nfolds
            trainIdx = cv.training(k); testIdx = cv.test(k);
            mdl = fitglm(SpkMat(trainIdx,:), logical(odorIDs(trainIdx)),'Distribution', 'binomial');
            % mdl=fitcdiscr(SpkMat(trainIdx,:), odorIDs(trainIdx));
            % predict regression output
            Y_hat = predict(mdl, SpkMat(testIdx,:));
            pred = round(Y_hat); %glm outputs very small values instead of zeros (binomial) for some reason, round to zero
            ids = odorIDs(testIdx);
            fract_correct(k) =nanmean(pred == ids);
        end
        fract_correct_real = mean(fract_correct);
        
        %Do shuffling
        nBoots=100; wb=waitbar(0,'Starting Boots');
        for boot = 1:nBoots
            waitbar(boot/nBoots,wb,'Running Boots')
            shuff_odorIDs=odorIDs(randperm(length(odorIDs))); % reorder the matrix so it doesnt match odors
            for k=1:nfolds
                trainIdx = cv.training(k); testIdx = cv.test(k);
                mdl = fitglm(SpkMat(trainIdx,:), shuff_odorIDs(trainIdx),'Distribution', 'binomial');
                % mdl=fitcdiscr(SpkMat(trainIdx,:), shuff_odorIDs(trainIdx));
                Y_hat_shuff = predict(mdl, SpkMat(testIdx,:));
                pred_shuff = round(Y_hat_shuff); %glm outputs very small values instead of zeros (binomial) for some reason, round to zero
                ids_shuff = shuff_odorIDs(testIdx);
                fract_correct_shuff(k) = nanmean(pred_shuff == ids_shuff);
            end
            fract_correct_null(boot) = mean(fract_correct_shuff);
        end
        
        close(wb);
        p_correct=1-normcdf(fract_correct_real,nanmean(fract_correct_null),...
            nanstd(fract_correct_null));
        fprintf('%s %d %s odor decoding trials:%d, correct%.f%%  null=%.f%% p=%.3f\n',...
            SuperRat(ses).name,SuperRat(ses).daynum,region{i},length(odorIDs),...
            fract_correct_real*100,nanmean(fract_correct_null)*100,p_correct);
        if verbose
            kill; figure;
            histogram(fract_correct_null,25); hold on; y=get(gca,'YLim');
            plot(repmat(fract_correct_real,1,2),[0,y(2)],'r');
            plot(repmat(nanmean(fract_correct_null),1,2),[0,y(2)],'b');
            title(sprintf('%s day %d %s',SuperRat(ses).name,SuperRat(ses).daynum, region{i}));
            legend('null','Odor Decode mean','null mean');
            drawnow;
        end
        
        decoderprobs(ses,i,1)=fract_correct_real;
        decoderprobs(ses,i,2)=nanmean(fract_correct_null);
    end
    % this method of taking the mean run probably doesnt make much sense
    % i think i'd rather do something like a matched timewindow analysis
    
end
warning('on','all');


%%
figure;
for i=1:2
    subplot(1,2,i);
    errorbar([1 2],[nanmean(decoderprobs(:,i,1)) nanmean(decoderprobs(:,i,2))],...
        [SEM(decoderprobs(:,i,1)) SEM(decoderprobs(:,i,1))],'.');
    xlim([.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'real','null'});
    title(sprintf('%s odor decoder',region{i}));
end
%}

%% this is the same decoder but now also decoding odor on track runs
% we also arent using claires decoder here

%{
%{

algorithm for decoding the runs:
the question:
    -where along the run does the neural ensemble resemble that of past
    odor coding?
       -maybe this approach is to pull the odor ensemble, then correlate it
       to the following timebins for each 800 msec run, and then aggregate
       across
    -does the odor discriminating code persist from odor coding into the
    delay period?
%}
veltimesmooth=8; speedthreshold=3;
verbose=0;
nBoots=200;
nextrunonly=1;
winsize=1;
mytemplate={sprintf('%d boots',nBoots),'Odor Decode','Run Decode';'CA1',[],[];'PFC',[],[]};
mytemplate(:,:,2)=mytemplate; mytemplate{1}='Success Prob real, null';
nfolds = 5;
warning('off','all');
savedir=uigetdir;

for ses=1:length(SuperRat)
    
    % take only correct trials during active epochs
    odorOK=ismember(SuperRat(ses).trialdata.EpochInds(:,2),SuperRat(ses).RunEpochs);
    odorflip=~logical(SuperRat(ses).trialdata.CorrIncorr10(odorOK)); % wrong trials
    odorStarts=SuperRat(ses).trialdata.sniffstart(odorOK); % start time
    odorIDs=SuperRat(ses).trialdata.leftright10(odorOK); % odor id
    badtrials=isnan(odorIDs);
    % somehow nans snuck into the odor ids? wtf?
    odorStarts(badtrials)=[]; odorIDs(badtrials)=[]; odorflip(badtrials)=[];
    odorIDs(odorflip)=double(~odorIDs(odorflip)); % flip the sign of the wrong trials
    %
    odorSpkMat=[]; foundOdorspikes={}; % grab big old matrix
    for j=1:length(SuperRat(ses).units)
        [~,odorSpkMat(:,j),~,~,~,foundOdorspikes(:,j)]=event_spikes(SuperRat(ses).units(j).ts,...
            odorStarts,0,winsize);
    end
    % remove cells that dont contribute to odor coding
    badunits=sum(odorSpkMat>0,1)<10; % claire actualy kills cells who are active for fewer than 10 trials...
    
    % now gather all the runs.
    fulltraj=SuperRat(ses).LinCoords;
    
    % first, get a true velocity, filter (1 if running, 0 if not)
    % fulltraj(:,10)=SmoothMat2(fulltraj(:,7),[0 50],veltimesmooth)>=speedthreshold;
    % this might be a good place to figure out which cells are silent
    % for which epochs
    trajinds=[3 1; 3 2]; % only outbound this time
    keepinds=ismember(fulltraj(:,4:5),trajinds,'rows');
    temptraj=sortrows(fulltraj(keepinds,:),1);
    % temptraj(temptraj(:,10)==0)=[]; % now cut out when he's slow
    breaks=find(diff(temptraj(:,1))>1); % only use time breaks that last more than a second
    % concatenate the start of each traj, its end, and the origin and
    % destination at that first index
    runepochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1);...
        temptraj(end,1)] temptraj([1; breaks+1],[4 5])];
    runinds=[[1; breaks+1] [breaks; size(temptraj,1)] temptraj([1; breaks+1],[4 5])];
    runepochs(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
    runinds(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
    
    % here would be the place to pick the runs that follow each correct
    % odor sample...
    if nextrunonly
        okruns=[];
        for rn=1:length(odorStarts)
            % find the next run start
            nextind=find(runepochs(:,1)>odorStarts(rn,1),1,'first');
            OKrunepochs(rn,:)=runepochs(nextind,:);
            OKruninds(rn,:)=runinds(nextind,:);
        end
    else
        OKrunepochs=runepochs;
        OKruninds=runinds;
    end
    % identify each run epoch by which type of run it is
    % now get the likelihood of each spike per msec
    runIDs=ismember(OKruninds(:,[3 4]),trajinds(1,:),'rows');
    runSpkMat=[]; foundRunSpikes={};
    for j=1:length(SuperRat(ses).units)
        [~,runSpkMat(:,j),~,~,~,foundRunSpikes{:,j}]=event_spikes(SuperRat(ses).units(j).ts,...
            OKrunepochs(:,1),0,OKrunepochs(:,2)-OKrunepochs(:,1));
    end
    
    % now run decoder
    thisdecode=mytemplate;
    if verbose, kill; figure; end
    
    for i=1:length(region)
        inRegion=cellfun(@(a) contains(a,region{i},'IgnoreCase',true), {SuperRat(ses).units.area});
        SpkMat=odorSpkMat(:,inRegion & ~badunits); % p(spike| 1 ms timebin)
        
        % true decoding
        cv = cvpartition(length(odorIDs), 'kfold',nfolds);
        for k=1:nfolds
            trainIdx = cv.training(k); testIdx = cv.test(k);
            mdl = fitglm(SpkMat(trainIdx,:), odorIDs(trainIdx),'Distribution', 'binomial');
            % predict regression output
            Y_hat = predict(mdl, SpkMat(testIdx,:));
            pred = round(Y_hat); %glm outputs very small values instead of zeros (binomial) for some reason, round to zero
            ids = odorIDs(testIdx);
            fract_correct(k) = sum(pred == ids)/length(pred);
        end
        fract_correct_real = mean(fract_correct);
        
        %Do shuffling
        wb=waitbar(0,'Starting Boots');
        for boot = 1:nBoots
            waitbar(boot/nBoots,wb,'Running Boots')
            shuff_odorIDs=odorIDs(randperm(length(odorIDs))); % reorder the matrix so it doesnt match odors
            for k=1:nfolds
                trainIdx = cv.training(k); testIdx = cv.test(k);
                mdl = fitglm(SpkMat(trainIdx,:), shuff_odorIDs(trainIdx),'Distribution', 'binomial');
                Y_hat_shuff = predict(mdl, SpkMat(testIdx,:));
                pred_shuff = round(Y_hat_shuff);
                ids_shuff = shuff_odorIDs(testIdx);
                fract_correct_shuff(k) = sum(pred_shuff == ids_shuff)/length(pred_shuff);
            end
            fract_correct_null(boot) = mean(fract_correct_shuff);
        end
        close(wb);
        
        p_correct=1-normcdf(fract_correct_real,nanmean(fract_correct_null),...
            nanstd(fract_correct_null));
        
        % now decode runs
        Y_hat= predict(mdl, runSpkMat(:,inRegion & ~badunits));
        pred = round(Y_hat); %glm outputs very small values instead of zeros (binomial) for some reason, round to zero
        fract_correct_run = sum(pred == runIDs)/length(pred);
        bino_p_run=1-binocdf(sum(pred==runIDs),length(pred),.5);
        % or use your null
        boot_p_run=1-normcdf(fract_correct_run,nanmean(fract_correct_null),...
            nanstd(fract_correct_null));
        if verbose
            subplot(1,2,i);
            histogram(fract_correct_null,25); hold on; y=get(gca,'YLim');
            od=plot(repmat(fract_correct_real,1,2),[0,y(2)],'r');
            rd=plot(repmat(fract_correct_run,1,2),[0,y(2)],'m');
            title(sprintf('%s day %d %s',SuperRat(ses).name,SuperRat(ses).daynum, region{i}));
            legend([od rd],'Odor Decode','Run Decode');
            drawnow;
        end
        
        fprintf('%s %d %s trials %d null: %.f%%, decode odor (%d%%) p=%.4f, run decode (%d%%) bino p=%.4f boot p=%.2e \n',...
            SuperRat(ses).name, SuperRat(ses).daynum, region{i}, length(odorIDs),...
            nanmean(fract_correct_null)*100,round(fract_correct_real*100),...
            p_correct, round(fract_correct_run*100), bino_p_run, boot_p_run);
        
        % ought to save it somewhere...
        thisdecode{i+1,2,1}=[fract_correct_real nanmean(fract_correct_null)];
        thisdecode{i+1,3,1}=[fract_correct_run nanmean(fract_correct_null)];
        thisdecode{i+1,2,2}=p_correct; thisdecode{i+1,3,2}=boot_p_run;
        
    end
    SuperRat(ses).DecodePvals=thisdecode;
end

if any(savedir~=0)
    today=datestr(now); mydate=today(1:find(today==' ',1,'first')-1);
    save(fullfile(savedir,sprintf('ClaireData%s',mydate)),'SuperRat');
    fprintf('Finished Boot and Saved out \n');
else
    fprintf('finished boot, didnt save \n');
end
%% lets just plot these performances.. i'm not sure they're very good
figure;
for i =1:2
    
    realProbs=cellfun(@(a) a{i+1,2,1}(1), {SuperRat.DecodePvals});
    nullProbs=cellfun(@(a) a{i+1,2,1}(2), {SuperRat.DecodePvals});
    
    subplot(2,2,i);
    errorbar([1 2],[nanmean(realProbs), nanmean(nullProbs)],[SEM(realProbs),SEM(nullProbs)],'.');
    xlim([.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'real','null'});
    title(sprintf('%s odor decoder',region{i}));
    
    realProbs=cellfun(@(a) a{i+1,3,1}(1), {SuperRat.DecodePvals});
    nullProbs=cellfun(@(a) a{i+1,3,1}(2), {SuperRat.DecodePvals});
    
    subplot(2,2,i+2);
    errorbar([1 2],[nanmean(realProbs), nanmean(nullProbs)],[SEM(realProbs),SEM(nullProbs)],'.');
    title(region{i}); xlim([.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'real','null'});
    title(sprintf('%s run decoder',region{i}));
    
end
linkaxes(get(gcf,'Children'),'y')

%}

%% now using the bayesian method with poisson firing
% this uses a training set of odor spikes, and decodes left out odor trials
% as well as run trials.  It may be a good idea to do the same leave one
% out for this though...

% quartiles or quintiles?
% four
runpos=[1 25; 25 50; 50 75; 75 100]; 
% or five
runpos=[1 20; 21 40; 41 60; 61 80; 81 100]; 

% a quick wrapper



% calculate the probability of each time given this activity
rstream = RandStream('mt19937ar','Seed',16);
RandStream.setGlobalStream(rstream);

%p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
onlyodorruns=1;
decoderprobs=[];
warning('off','all');
winsize=1; % second of odor sampling
% this is about 0.25 seconds of buffer time between end sample and start
% run
deadspace=0.5; % seconds following odor sampling start to never consider
verbose=0;
nDolds = 5;
nBoots=200;
odorNames=[1 0];
veltimesmooth=8;
speedthreshold=3;
region={'CA1','PFC'};

mytemplate={sprintf('%d boots',nBoots),'Odor Decode','Run Decode';'CA1',[],[];'PFC',[],[]};
mytemplate(:,:,2)=mytemplate; mytemplate{1}='Success Prob real, null';
for ses=1:length(SuperRat)
    % take only correct trials during active epochs
    odorOK=ismember(SuperRat(ses).trialdata.EpochInds(:,2),SuperRat(ses).RunEpochs);
    odorflip=~logical(SuperRat(ses).trialdata.CorrIncorr10(odorOK)); % wrong trials
    odorStarts=SuperRat(ses).trialdata.sniffstart(odorOK); % start time
    odorEnds=SuperRat(ses).trialdata.sniffend(odorOK); % start time
    odorIDs=SuperRat(ses).trialdata.leftright10(odorOK); % odor id
    badtrials=isnan(odorIDs);
    % somehow nans snuck into the odor ids? wtf?
    odorStarts(badtrials)=[]; odorEnds(badtrials)=[];
    odorIDs(badtrials)=[]; odorflip(badtrials)=[];
    odorIDs(odorflip)=double(~odorIDs(odorflip)); % flip the sign of the wrong trials
    %
    odorSpkMat=[]; foundOdorspikes={}; % grab big old matrix
    for j=1:length(SuperRat(ses).units)
        [~,odorSpkMat(:,j),~,~,~,foundOdorspikes(:,j)]=event_spikes(SuperRat(ses).units(j).ts,...
            odorStarts,odorEnds-odorStarts,winsize);
    end
    % remove cells that dont contribute to odor coding
    badunits=sum(odorSpkMat>0,1)<10 | nanmean(odorSpkMat,1)>9; % claire actualy kills cells who are active for fewer than 10 trials...
    
    % now gather all the runs.
    fulltraj=SuperRat(ses).LinCoords;
    
    % first, get a true velocity, filter (1 if running, 0 if not)
    fulltraj(:,10)=SmoothMat2(fulltraj(:,7),[0 50],veltimesmooth)>=speedthreshold;
    % this might be a good place to figure out which cells are silent
    % for which epochs
    trajinds=[3 1; 3 2]; % only outbound this time
    keepinds=ismember(fulltraj(:,4:5),trajinds,'rows');
    fulltraj=sortrows(fulltraj(keepinds,:),1);
    
    %
    % work up the tracking data into track bins
    %
    runSpkMat={}; runIDs={};
    for ile=1:length(runpos)
        % only use times when the linear trajectory is between...
        okpos=runpos(ile,:);% segment out the arm
        oktrackpos=fulltraj(:,8)>okpos(1) & fulltraj(:,8)<okpos(2);
        temptraj=fulltraj(oktrackpos,:);
        
        killslow=0; % kill slow moving times
        if killslow
            temptraj(temptraj(:,10)==0)=[]; % now cut out when he's slow
        end
        
        % kill bins that are too close to the sampling period (1 sec)
        
        if deadspace>0
            for i=1:length(odorEnds)
                % kill the n seconds after the odor end
                timediffs=temptraj(:,1)-odorEnds(i,1);
                % for that end, kill all following tracking coords until
                % deadspace ends
                badtimes=timediffs<deadspace & timediffs>0;
                temptraj(badtimes,:)=[];
                hadtocut(i)=sum(badtimes);
            end
        end
        
       
        breaks=find(diff(temptraj(:,1))>1); % only use time breaks that last more than a second
        % concatenate the start of each traj, its end, and the origin and
        % destination at that first index
        runepochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1);...
            temptraj(end,1)] temptraj([1; breaks+1],[4 5])];
        runinds=[[1; breaks+1] [breaks; size(temptraj,1)] temptraj([1; breaks+1],[4 5])];
        %runepochs(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        %runinds(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        
        % here would be the place to pick the runs that follow each correct
        % odor sample...
        
        OKruninds=[]; OKrunepochs=[];
        if onlyodorruns
            for rn=1:length(odorStarts)
                % find the next run start
                nextind=find(runepochs(:,1)>odorStarts(rn,1),1,'first');
                OKrunepochs(rn,:)=runepochs(nextind,:);
                OKruninds(rn,:)=runinds(nextind,:);
            end
        else
            OKrunepochs=runepochs;
            OKruninds=runinds;
        end
        
        % identify each run epoch by which type of run it is
        % now get the likelihood of each spike per msec
        runIDs{ile}=ismember(OKruninds(:,[3 4]),trajinds(1,:),'rows');

        for j=1:length(SuperRat(ses).units)
            [~,runSpkMat{ile}(:,j)]=event_spikes(SuperRat(ses).units(j).ts,...
                OKrunepochs(:,1),0,OKrunepochs(:,2)-OKrunepochs(:,1));
        end
    end


    thisdecode=mytemplate;
    for i=1:length(region)
        inRegion=cellfun(@(a) contains(a,region{i},'IgnoreCase',true), {SuperRat(ses).units.area});
        % here is where you can use ONLY PLACE CELLS
        if sum(inRegion)>5
            SpkMat=odorSpkMat(:,inRegion & ~badunits)*winsize; % p(spike| 1 ms timebin)
            
            % run the within group decoder
            [p_correct,fract_correct_real,fract_correct_null]...
                = nFoldBayesPoisson(SpkMat,odorIDs,5,200);

            thisdecode{i+1,2,1}=[fract_correct_real nanmean(fract_correct_null)];
            thisdecode{i+1,2,2}=p_correct;
            
            % now decode runs based on odor period
            odorPriors=mean(SpkMat(odorIDs==1,:)); % 1,nunits,1
            odorPriors(:,:,2)=mean(SpkMat(odorIDs==0,:)); % 1,nunits,2
            for ile=1:length(runpos)
                testmat=runSpkMat{ile}(:,inRegion & ~badunits); % ntrials x nunits
                trueOdors=runIDs{ile}; % true run sides
                prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
                summat=sum(odorPriors,2); % sum across units (mean_i)
                tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
                probmat=prodmat.*exp(-tau.*summat);
                % sum p across potential odors must ==1 (the c term)
                realProb=probmat./sum(probmat,3);
                [~,decoded]=max(squeeze(realProb),[],2);
                fract_correct_run=nanmean(trueOdors==odorNames(decoded)');
                
                            
                for boot=1:200
                    scrambledecode=decoded(randperm(length(decoded)));
                    fract_correct_shuff(boot)=nanmean(trueOdors==odorNames(scrambledecode)');
                end
                fract_correct_null=nanmean(fract_correct_shuff);
                
                p_run=1-normcdf(fract_correct_run,nanmean(fract_correct_shuff),...
                    nanstd(fract_correct_shuff));

                thisdecode{i+1,3,1}=[thisdecode{i+1,3,1}; [fract_correct_run nanmean(fract_correct_null)]];
                thisdecode{i+1,3,2}(ile)=p_run;
            end

            fprintf('%s %d %s trials:%d,  odor decoding correct%.f%%, run decoding correct%.f%%,  null=%.f%% p=%.3f, p%.3f\n',...
                SuperRat(ses).name,SuperRat(ses).daynum,region{i},length(odorIDs),...
                fract_correct_real*100, fract_correct_run*100,...
                nanmean(fract_correct_null)*100,p_correct,p_run);
            
        else

            thisdecode{i+1,2,1}=[nan nan];
            thisdecode{i+1,3,1}=[nan nan];
            thisdecode{i+1,2,2}=nan; thisdecode{i+1,3,2}=nan;
            fprintf('%s %d %s not used, not enough units... \n',...
                SuperRat(ses).name,SuperRat(ses).daynum,region{i});
        end
    end
    SuperRat(ses).BayesDecodePvals=thisdecode;
    % this method of taking the mean run probably doesnt make much sense
    % i think i'd rather do something like a matched timewindow analysis
    
end
warning('on','all');



%% work up data for odor period to odor period

figure;
for i=1:2
    
    realProbs=cellfun(@(a) a{i+1,2,1}(1), {SuperRat.BayesDecodePvals});
    nullProbs=cellfun(@(a) a{i+1,2,1}(2), {SuperRat.BayesDecodePvals});
    
    subplot(2,2,i);
    errorbar([1 2],[nanmean(realProbs), nanmean(nullProbs)],[nanstd(realProbs),nanstd(nullProbs)],...
        '.','CapSize',0);
    xlim([.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'real','null'});
    title(sprintf('%s odor decoder',region{i}));
    % could get away with a tiedrank test here on the diffs
    [a,b]=ranksum(realProbs,nullProbs);

    fprintf('%s odor decoder P= %.3f \n',region{i},a); 
    
    realProbs=cellfun(@(a) nanmean(a{i+1,3,1}(:,1)), {SuperRat.BayesDecodePvals});
    nullProbs=cellfun(@(a) nanmean(a{i+1,3,1}(:,2)), {SuperRat.BayesDecodePvals});
    
    subplot(2,2,i+2);
    errorbar([1 2],[nanmean(realProbs), nanmean(nullProbs)],[SEM(realProbs),SEM(nullProbs)],...
        '.','CapSize',0);
    title(region{i}); xlim([.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'real','null'});
    title(sprintf('%s run decoder',region{i}));
    [a,b]=ranksum(realProbs,nullProbs);
    fprintf('%s run decoder P= %.3f \n',region{i},a);
end


linkaxes(get(gcf,'Children'),'y')
%% work up data for odor period to track periods
figure;
linecolors=[lines(2); [.7 .7 .7]];
for i=1:2
    realProbs=cellfun(@(a) a{i+1,3,1}(:,1), {SuperRat.BayesDecodePvals}, 'UniformOutput', false);
    realProbMat=cell2mat(realProbs(cellfun(@(a) length(a)>1,realProbs)));
    nullProbs=cellfun(@(a) a{i+1,3,1}(:,2), {SuperRat.BayesDecodePvals},'UniformOutput', false);
    nullProbMat=cell2mat(nullProbs(cellfun(@(a) length(a)>1,nullProbs)));

    subplot(1,2,i);
    for ile=1:size(realProbMat,1)
        ho=errorbar(ile,nanmean(realProbMat(ile,:)),nanstd(realProbMat(ile,:)),...
            '.','Color',linecolors(i,:),'CapSize',0);
        hold on;
        ha=errorbar(ile+.1,nanmean(nullProbMat(ile,:)),nanstd(nullProbMat(ile,:)),...
            '.','Color',linecolors(3,:),'CapSize',0);
        [p,h]=signrank(realProbMat(ile,:)-nullProbMat(ile,:));
        %[p,h] = ranksum(realProbMat(ile,:), nullProbMat(ile,:));
        if p<.05
           text(ile,max(nanmean(realProbMat,2))*1.1,'*','FontSize',16);
        end
        fprintf('%s quintile %d p=%.4f \n',region{i},ile,p);
           
    end
    legend('Real Data','Null');
    xlim([.5 5.5]);
    xlabel('track quintile'); ylabel('Decoding success probability');
end
linkaxes;
sgtitle('Odor coding ensembles persist into track run');

    
%% and now a control, how easy is it to decode the run after training on the run?
% e.g. train on a segment and test on THAT SEGMENT

% four
runpos=[1 25; 25 50; 50 75; 75 100]; 
% or five
runpos=[1 20; 21 40; 41 60; 61 80; 81 100]; 



% calculate the probability of each time given this activity
rstream = RandStream('mt19937ar','Seed',16);
RandStream.setGlobalStream(rstream);

%p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
onlyodorruns=1;
decoderprobs=[];
warning('off','all');
winsize=1; % second of odor sampling
% this is about 0.25 seconds of buffer time between end sample and start
% run
deadspace=0.5; % seconds following odor sampling start to never consider
verbose=0;
nDolds = 5;
nBoots=200;
odorNames=[1 0];
veltimesmooth=8;
speedthreshold=3;
region={'CA1','PFC'};

mytemplate={sprintf('%d boots',nBoots),'Run-Run Decode';'CA1',[];'PFC',[]};
mytemplate(:,:,2)=mytemplate; mytemplate{1}='Success Prob real, null';

for ses=1:length(SuperRat)
   % take only correct trials during active epochs
    odorOK=ismember(SuperRat(ses).trialdata.EpochInds(:,2),SuperRat(ses).RunEpochs);
    odorStarts=SuperRat(ses).trialdata.sniffstart(odorOK); % start time
    odorEnds=SuperRat(ses).trialdata.sniffend(odorOK); % End time
    odorflip=~logical(SuperRat(ses).trialdata.CorrIncorr10(odorOK)); % wrong trials

    odorIDs=SuperRat(ses).trialdata.leftright10(odorOK); % odor id
    badtrials=isnan(odorIDs);
    % somehow nans snuck into the odor ids? wtf?
    odorStarts(badtrials)=[]; odorEnds(badtrials)=[];
    odorIDs(badtrials)=[]; odorflip(badtrials)=[];
    odorIDs(odorflip)=double(~odorIDs(odorflip)); % flip the sign of the wrong trials
    % now gather all the runs.
    fulltraj=SuperRat(ses).LinCoords;
    
    % first, get a true velocity, filter (1 if running, 0 if not)
    fulltraj(:,10)=SmoothMat2(fulltraj(:,7),[0 50],veltimesmooth)>=speedthreshold;
    % this might be a good place to figure out which cells are silent
    % for which epochs
    trajinds=[3 1; 3 2]; % only outbound this time
    keepinds=ismember(fulltraj(:,4:5),trajinds,'rows');
    fulltraj=sortrows(fulltraj(keepinds,:),1);
    
    %
    % work up the tracking data into track bins
    %
    runSpkMat={}; runIDs={};
    for ile=1:length(runpos)
        % only use times when the linear trajectory is between...
        okpos=runpos(ile,:);% segment out the arm
        oktrackpos=fulltraj(:,8)>okpos(1) & fulltraj(:,8)<okpos(2);
        temptraj=fulltraj(oktrackpos,:);
        
        killslow=0; % kill slow moving times
        if killslow
            temptraj(temptraj(:,10)==0)=[]; % now cut out when he's slow
        end
        
        % kill bins that are too close to the sampling period (1 sec)
        if deadspace>0
            for i=1:length(odorEnds)
                % Seconds after each end
                timediffs=temptraj(:,1)-odorEnds(i,1);
                afterstart=temptraj(:,1)-odorStarts(i,1);
                % add dead space, kill anything after the end but within
                % dead space time
                badtimes=timediffs<deadspace & afterstart>0;
                temptraj(badtimes,:)=[];
                hadtocut(i)=sum(badtimes);
            end
        end
        
        % make this 50 to remove distal arm timestamps
        killsidearm=0;
        if killsidearm>0
            temptraj(temptraj(:,8)>killsidearm,:)=[];
        end
        
        
        breaks=find(diff(temptraj(:,1))>1); % only use time breaks that last more than a second
        % concatenate the start of each traj, its end, and the origin and
        % destination at that first index
        runepochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1);...
            temptraj(end,1)] temptraj([1; breaks+1],[4 5])];
        runinds=[[1; breaks+1] [breaks; size(temptraj,1)] temptraj([1; breaks+1],[4 5])];
        %runepochs(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        %runinds(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        
        % here would be the place to pick the runs that follow each correct
        % odor sample...
        
        OKruninds=[]; OKrunepochs=[];
        if onlyodorruns
            for rn=1:length(odorStarts)
                % find the next run start
                nextind=find(runepochs(:,1)>odorStarts(rn,1),1,'first');
                OKrunepochs(rn,:)=runepochs(nextind,:);
                OKruninds(rn,:)=runinds(nextind,:);
            end
        else
            OKrunepochs=runepochs;
            OKruninds=runinds;
        end
        
        % identify each run epoch by which type of run it is
        % now get the likelihood of each spike per msec
        runIDs{ile}=ismember(OKruninds(:,[3 4]),trajinds(1,:),'rows');

        for j=1:length(SuperRat(ses).units)
            [~,runSpkMat{ile}(:,j)]=event_spikes(SuperRat(ses).units(j).ts,...
                OKrunepochs(:,1),0,OKrunepochs(:,2)-OKrunepochs(:,1));
        end
    end


    thisdecode=mytemplate;
    for i=1:length(region)
        inRegion=cellfun(@(a) contains(a,region{i},'IgnoreCase',true), {SuperRat(ses).units.area});
        if sum(inRegion)>5
            for ile=1:length(runpos)
                tempMat=runSpkMat{ile}(:,inRegion); % ntrials x nunits
                badunits=sum(tempMat>0)<10; % spikes on fewer than 10 runs
                tempMat(:,badunits)=[];
                trueOdors=runIDs{ile}; % true run sides
                
                % run the within group decoder
                [p_correct,fract_correct_real,fract_correct_null]...
                    = nFoldBayesPoisson(tempMat,trueOdors,5,200);
                
                thisdecode{i+1,2,1}(ile,:)=[fract_correct_real nanmean(fract_correct_null)];
                thisdecode{i+1,2,2}(ile)=p_correct;
                
                
                
                fprintf('%s %d %s routepiece %d, trials:%d, run decoding correct%.f%%,  null=%.f%% p=%.3f,\n',...
                    SuperRat(ses).name,SuperRat(ses).daynum,region{i},ile,length(trueOdors),...
                    fract_correct_real*100,nanmean(fract_correct_null)*100,p_correct);
            end
        else

            thisdecode{i+1,2,1}=[nan nan];
            thisdecode{i+1,3,1}=[nan nan];
            thisdecode{i+1,2,2}=nan; thisdecode{i+1,3,2}=nan;
            fprintf('%s %d %s not used, not enough units... \n',...
                SuperRat(ses).name,SuperRat(ses).daynum,region{i});
        end
    end
    SuperRat(ses).BayesDecodePvals_run_run=thisdecode;
    % this method of taking the mean run probably doesnt make much sense
    % i think i'd rather do something like a matched timewindow analysis
    
end
warning('on','all');
%% and analyze the results...

figure;
linecolors=[lines(2); [.7 .7 .7]];
for i=1:2
    realProbs=cellfun(@(a) a{i+1,2,1}(:,1), {SuperRat.BayesDecodePvals_run_run}, 'UniformOutput', false);
    realProbMat=cell2mat(realProbs(cellfun(@(a) length(a)>1,realProbs)));
    nullProbs=cellfun(@(a) a{i+1,2,1}(:,2), {SuperRat.BayesDecodePvals_run_run},'UniformOutput', false);
    nullProbMat=cell2mat(nullProbs(cellfun(@(a) length(a)>1,nullProbs)));

    subplot(1,2,i);
    for ile=1:size(realProbMat,1)
        ho=errorbar(ile-.025,nanmean(realProbMat(ile,:)),nanstd(realProbMat(ile,:)),...
            '.','Color',linecolors(i,:),'CapSize',0);
        hold on;
        ha=errorbar(ile+.025,nanmean(nullProbMat(ile,:)),nanstd(nullProbMat(ile,:)),...
            '.','Color',linecolors(3,:),'CapSize',0);
        [p,h]=signrank(realProbMat(ile,:)-nullProbMat(ile,:));
        %[p,h] = ranksum(realProbMat(ile,:), nullProbMat(ile,:));
        if p<.05
            if p<0.01
                text(ile-.25,max(nanmean(realProbMat,2))*1.1,'**','FontSize',16);
            else
                text(ile-.1,max(nanmean(realProbMat,2))*1.1,'*','FontSize',16);
            end
        end
    end
    legend('Real Data','Null');
    xlim([.5 5.5]); box off; ylim([.4 1]);
    xlabel('track quintile'); ylabel('Decoding success probability');
end
sgtitle('splitters are abundant across the track')
linkaxes;

%% now for each spatial bin test the performance on each other bin


% four
runpos=[1 25; 25 50; 50 75; 75 100]; 
% or five
runpos=[1 20; 21 40; 41 60; 61 80; 81 100]; 

%or eight
%runpos=[1 12; 13 24; 25 36; 37 50; 51 62; 63 74; 75 87; 88 100]; 
% a quick wrapper



% calculate the probability of each time given this activity
rstream = RandStream('mt19937ar','Seed',16);
RandStream.setGlobalStream(rstream);

%p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
onlyodorruns=1;
decoderprobs=[];
warning('off','all');
winsize=1; % second of odor sampling
% this is generally about 0.25 seconds of buffer time between end sample and start
% run
daytimer=tic;
deadspace=0.5; % seconds following odor sampling start to never consider
verbose=0;
nDolds = 5; % n folds for decoder
nBoots=200; % n boots for the null distr
odorNames=[1 0];
veltimesmooth=8; % seconds for the velocity smoothing kernel
speedthreshold=3; % lowest movement speed
region={'CA1','PFC'};

mytemplate={sprintf('null, %d boots',nBoots),'Run-Run Decode (study x test)';'CA1',[];'PFC',[]};
mytemplate(:,:,2)=mytemplate; mytemplate{1}='Success Prob real';

for ses=1:length(SuperRat)
    sesstimer=tic;
   % take only correct trials during active epochs
    odorOK=ismember(SuperRat(ses).trialdata.EpochInds(:,2),SuperRat(ses).RunEpochs);
    odorStarts=SuperRat(ses).trialdata.sniffstart(odorOK); % start time
    odorEnds=SuperRat(ses).trialdata.sniffend(odorOK); % End time
    odorflip=~logical(SuperRat(ses).trialdata.CorrIncorr10(odorOK)); % wrong trials

    odorIDs=SuperRat(ses).trialdata.leftright10(odorOK); % odor id
    badtrials=isnan(odorIDs);
    % somehow nans snuck into the odor ids? wtf?
    odorStarts(badtrials)=[]; odorEnds(badtrials)=[];
    odorIDs(badtrials)=[]; odorflip(badtrials)=[];
    odorIDs(odorflip)=double(~odorIDs(odorflip)); % flip the sign of the wrong trials
    % now gather all the runs.
    fulltraj=SuperRat(ses).LinCoords;
    
    % first, get a true velocity, filter (1 if running, 0 if not)
    fulltraj(:,10)=SmoothMat2(fulltraj(:,7),[0 50],veltimesmooth)>=speedthreshold;
    % this might be a good place to figure out which cells are silent
    % for which epochs
    trajinds=[3 1; 3 2]; % only outbound this time
    keepinds=ismember(fulltraj(:,4:5),trajinds,'rows');
    fulltraj=sortrows(fulltraj(keepinds,:),1);
    
    %
    % work up the tracking data into track bins
    %
    runSpkMat={}; runIDs={};
    for ile=1:length(runpos)
        % only use times when the linear trajectory is between...
        okpos=runpos(ile,:);% segment out the arm
        oktrackpos=fulltraj(:,8)>okpos(1) & fulltraj(:,8)<okpos(2);
        temptraj=fulltraj(oktrackpos,:);
        
        killslow=0; % kill slow moving times
        if killslow
            temptraj(temptraj(:,10)==0)=[]; % now cut out when he's slow
        end
        
        % kill bins that are too close to the sampling period (1 sec)
        if deadspace>0
            for i=1:length(odorEnds)
                % Seconds after each end
                timediffs=temptraj(:,1)-odorEnds(i,1);
                afterstart=temptraj(:,1)-odorStarts(i,1);
                % add dead space, kill anything after the end but within
                % dead space time
                badtimes=timediffs<deadspace & afterstart>0;
                temptraj(badtimes,:)=[];
                hadtocut(i)=sum(badtimes);
            end
        end
        
        % make this 50 to remove distal arm timestamps
        killsidearm=0;
        if killsidearm>0
            temptraj(temptraj(:,8)>killsidearm,:)=[];
        end
        
        
        breaks=find(diff(temptraj(:,1))>1); % only use time breaks that last more than a second
        % concatenate the start of each traj, its end, and the origin and
        % destination at that first index
        runepochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1);...
            temptraj(end,1)] temptraj([1; breaks+1],[4 5])];
        runinds=[[1; breaks+1] [breaks; size(temptraj,1)] temptraj([1; breaks+1],[4 5])];
        %runepochs(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        %runinds(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        
        % here would be the place to pick the runs that follow each correct
        % odor sample...
        
        OKruninds=[]; OKrunepochs=[];
        if onlyodorruns
            for rn=1:length(odorStarts)
                % find the next run start
                nextind=find(runepochs(:,1)>odorStarts(rn,1),1,'first');
                if ~isempty(nextind) % sometimes the last odor doesnt have a following run
                    OKrunepochs(rn,:)=runepochs(nextind,:);
                    OKruninds(rn,:)=runinds(nextind,:);
                end
            end
        else
            OKrunepochs=runepochs;
            OKruninds=runinds;
        end
        
        % identify each run epoch by which type of run it is
        % now get the likelihood of each spike per msec
        runIDs{ile}=ismember(OKruninds(:,[3 4]),trajinds(1,:),'rows');

        for j=1:length(SuperRat(ses).units)
            [~,runSpkMat{ile}(:,j)]=event_spikes(SuperRat(ses).units(j).ts,...
                OKrunepochs(:,1),0,OKrunepochs(:,2)-OKrunepochs(:,1));
        end
        
        
    end
    
    % reshape all matrices to be the same size
    matsz=cellfun(@(a) size(a,1), runSpkMat);
    runSpkMat=cellfun(@(a) a(1:min(matsz),:), runSpkMat, 'UniformOutput', false);
    runIDs=cellfun(@(a) a(1:min(matsz)), runIDs, 'UniformOutput', false);
    
    %
    % run decoder
    %

    thisdecode=mytemplate;
    for i=1:length(region)
        
        %%%%%%% FILTER UNITS HERE %%%%%%
        % right now its a place field and in region filter
        inRegion=cellfun(@(a) contains(a,region{i},'IgnoreCase',true), {SuperRat(ses).units.area}) & ...       
        cell2mat(cellfun(@(a) any(a), {SuperRat(ses).units.FiresDuringRun}, 'UniformOutput',false));
        %cell2mat(cellfun(@(a) any(a(1:2)), {SuperRat(ses).units.PFexist},'uniformOutput',false));
        
        if sum(inRegion)>5
            for ile=1:length(runpos)
                for ile2=1:length(runpos)
                tempMat=runSpkMat{ile}(:,inRegion); % ntrials x nunits
                tempMat2=runSpkMat{ile2}(:,inRegion);
                badunits=sum(tempMat>0)<10 | sum(tempMat2>0)<10; % spikes on fewer than 10 runs
                tempMat(:,badunits)=[]; tempMat2(:,badunits)=[];
                trueOdors=runIDs{ile}; % true run sides
                
                % run the within group decoder
                [p_correct,fract_correct_real,fract_correct_null]...
                    = nFoldBayesPoisson2(tempMat,tempMat2,trueOdors,nDolds,nBoots);
                
                thisdecode{i+1,2,1}(ile,ile2,:)=[fract_correct_real nanmean(fract_correct_null)];
                thisdecode{i+1,2,2}(ile,ile2)=p_correct;
                end
            end
        else
            thisdecode{i+1,2,1}=nan(length(runpos),length(runpos),2);
            thisdecode{i+1,2,2}=nan(length(runpos),length(runpos)); 
            fprintf('%s %d %s not used, not enough units... \n',...
                SuperRat(ses).name,SuperRat(ses).daynum,region{i});
        end
    end
    SuperRat(ses).BayesDecodePvals_run_run=thisdecode;
    % this method of taking the mean run probably doesnt make much sense
    % i think i'd rather do something like a matched timewindow analysis
    fprintf('Sess %d took %.2f minutes, overall we''re at %.2f minutes \n',...
        ses,toc(sesstimer)/60,toc(daytimer)/60)
end
warning('on','all');

%% and analyze the results
% this will now be in a 2d matrix, so we'll have to get clever about
% visualisations

figure; clear sp realProbMat
% hot can be part of autumn, and cold (pfc) can be mid of parula
for i=1:2
    % this is the fraction correct real (:,:,2) is p correct real
    realProbs=cellfun(@(a) permute(a{i+1,2,1}(:,:,1),[1 3 2]), {SuperRat.BayesDecodePvals_run_run}, 'UniformOutput', false);
    % so the shape is study set x session x test set
    realProbMat{i}=cell2mat(realProbs(cellfun(@(a) size(a,1)>1,realProbs)));
    nullProbs=cellfun(@(a) a{i+1,2,1}(:,:,2), {SuperRat.BayesDecodePvals_run_run}, 'UniformOutput', false);
    nullProbMat=cell2mat(cellfun(@(a) permute(a,[1 3 2]), nullProbs(cellfun(@(a) length(a)>1,nullProbs)),'UniformOutput', false));

    linestart=jet(size(realProbs{1},1)+6);
    linestart2=autumn(size(realProbs{1},1)+1);
    linecolors=[linestart2(1:end-1,:) ; [.7 .7 .7]];
    linecolors(:,:,2)=[linestart(2:end-5,:); [.7 .7 .7]];
    
    %linecolors=[parula(length(runpos)); [.7 .7 .7]];
    sp(i)=subplot(2,2,i);
    for ile=1:size(realProbMat{i},1)
        ho=errorbar((1:size(realProbMat{i},1))+ile*.05,nanmean(squeeze(realProbMat{i}(ile,:,:))),...
            nanstd(squeeze(realProbMat{i}(ile,:,:)),1),'.-','Color',linecolors(ile,:,i),'CapSize',0);
        hold on;
        ha=errorbar((1:size(realProbMat{i},1))-.02+ile*.05,nanmean(squeeze(nullProbMat(ile,:,:))),...
            nanstd(squeeze(nullProbMat(ile,:,:)),1),'.-','Color',linecolors(end,:,i),'CapSize',0);
        %[p,h]=signrank(realProbMat(ile,:)-nullProbMat(ile,:));
        %[p,h] = ranksum(realProbMat(ile,:), nullProbMat(ile,:));
       %if p<.05
            %if p<0.01
            %    text(ile,max(nanmean(realProbMat,2))*1.1,'**','FontSize',16);
           % else
            %    text(ile,max(nanmean(realProbMat,2))*1.1,'*','FontSize',16);
           % end
        %end
    end
    legend([ho ha],'Real Data','Null'); title(region{i});
    xlim([.5 5.5]); box off; ylim([.45 1]);
    xlabel('track quintile'); ylabel('Decoding success probability');
    myticks=get(gca,'XTickLabel');
    for k=1:length(myticks)
        myticks{k}=[' \color[rgb]{' num2str(squeeze(linecolors(k,:,i))) '} ' myticks{k}];
    end
    set(gca,'XTickLabel',myticks);
    
    
    
    %subplot(2,2,i+2);
    %imagesc(squeeze(nanmean(realProbMat,2)));
end
sgtitle('splitters are abundant across the track')
linkaxes(sp);

% statistics: friedman tests?
% so the question is 1. are there any differences in decoder success across
% training quintiles?
for i=1:2
    % first linearize the rows
    testmat=[];
    for tr=1:size(realProbMat{i},1)
        testmat(:,tr)=linearize(realProbMat{i}(tr,:,:));
    end
    testmat(sum(isnan(testmat),2)>0,:)=[];
    fprintf('\n %s any trainign quintile better than others? \n',region{i}); 
    [a,b,c]=friedman(testmat,size(realProbMat{i},1))
    if a<.01, d = multcompare(c), end
    
end


% and second, does the decoder do better when trained on the same quintile
% as the test set?
for i=1:2
    % first linearize the rows
    testmat={[],[]};
    for tr=1:size(realProbMat{i},1)
        testmat{1}=[testmat{1}; linearize(realProbMat{i}(tr,:,tr))];
        testmat{2}=[testmat{2}; linearize(realProbMat{i}(tr,:,(1:size(realProbMat{i},1))~=tr))];
    end

    fprintf('\n %s same vs other quintile \n',region{i});
    fprintf('n=%d \n',length(testmat{1}+length(testmat{2})))
    [a,b,c]=ranksum(testmat{1},testmat{2})
    %if a<.01, d = multcompare(c), end
    
end



%% an alternative visualisation of the data is....
% a matrix... but it looks worse...



%%     % so run the decoeer analysis on a run by run basis, so probably start each run

%{
odor selective neurons- of those cells, do they also do that during the
delay period - try to collapse along time for each run, and use that, dig
into when you can grab the


i think the question is about the route- are both regions splitting along
that route?

Maybe an answer is that on the route, we can decode position and route for
each region, and then we can test the relationship betwen odor/position
decoding and region decoding accuracy- if odor gets better decodes when
position does they're related, andif PFC decodes corerlate with CA1 decodes
thats interesting too...
     
%}

%% now using the  bayesian method with poisson firing
%{
% just a copy of the above
clear realProbs nullProbs;

runpos=[1 25; 25 50; 50 75; 75 100]; 

% a quick wrapper

for ile=1:length(runpos)

% calculate the probability of each time given this activity
rstream = RandStream('mt19937ar','Seed',16);
RandStream.setGlobalStream(rstream);

%p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
onlyodorruns=1;
decoderprobs=[];
warning('off','all');
winsize=1; % second of odor sampling
% this is about 0.25 seconds of buffer time between end sample and start
% run
deadspace=1.25; % seconds following odor sampling start to never consider
verbose=0;
nfolds = 5;
nBoots=200;
odorNames=[1 0];
veltimesmooth=8;
speedthreshold=3;
region={'CA1','PFC'};

mytemplate={sprintf('%d boots',nBoots),'Odor Decode','Run Decode';'CA1',[],[];'PFC',[],[]};
mytemplate(:,:,2)=mytemplate; mytemplate{1}='Success Prob real, null';
for ses=1:length(SuperRat)
   % take only correct trials during active epochs
    odorOK=ismember(SuperRat(ses).trialdata.EpochInds(:,2),SuperRat(ses).RunEpochs);
    odorflip=~logical(SuperRat(ses).trialdata.CorrIncorr10(odorOK)); % wrong trials
    odorStarts=SuperRat(ses).trialdata.sniffstart(odorOK); % start time
    odorIDs=SuperRat(ses).trialdata.leftright10(odorOK); % odor id
    badtrials=isnan(odorIDs);
    % somehow nans snuck into the odor ids? wtf?
    odorStarts(badtrials)=[]; odorIDs(badtrials)=[]; odorflip(badtrials)=[];
    odorIDs(odorflip)=double(~odorIDs(odorflip)); % flip the sign of the wrong trials
    %
    odorSpkMat=[]; foundOdorspikes={}; % grab big old matrix
    for j=1:length(SuperRat(ses).units)
        [~,odorSpkMat(:,j),~,~,~,foundOdorspikes(:,j)]=event_spikes(SuperRat(ses).units(j).ts,...
            odorStarts,0,winsize);
    end
    % remove cells that dont contribute to odor coding
    badunits=sum(odorSpkMat>0,1)<10 | nanmean(odorSpkMat,1)>9; % claire actualy kills cells who are active for fewer than 10 trials...
    
    % now gather all the runs.
    fulltraj=SuperRat(ses).LinCoords;
    
    % first, get a true velocity, filter (1 if running, 0 if not)
    fulltraj(:,10)=SmoothMat2(fulltraj(:,7),[0 50],veltimesmooth)>=speedthreshold;
    % this might be a good place to figure out which cells are silent
    % for which epochs
    trajinds=[3 1; 3 2]; % only outbound this time
    keepinds=ismember(fulltraj(:,4:5),trajinds,'rows');
    temptraj=sortrows(fulltraj(keepinds,:),1);
    
    %
    % work up the tracking data into track bins
    %
    runspikemat={}; runIDs={};

        % only use times when the linear trajectory is between...
        okpos=runpos(ile,:);% first segment of arm
        oktrackpos=temptraj(:,8)>okpos(1) & temptraj(:,8)<okpos(2);
        temptraj(~oktrackpos,:)=[];
        
        killslow=0; % kill slow moving times
        if killslow
            temptraj(temptraj(:,10)==0)=[]; % now cut out when he's slow
        end
        
        % kill bins that are too close to the sampling period (1 sec)
        
        if deadspace>0
            for i=1:length(odorStarts)
                % kill the n seconds after the odor end
                timediffs=temptraj(:,1)-odorStarts(i,1);
                badtimes=timediffs<deadspace & timediffs>0;
                temptraj(badtimes,:)=[];
                hadtocut(i)=sum(badtimes);
            end
        end
        
        % make this 50 to remove distal arm timestamps
        killsidearm=0;
        if killsidearm>0
            temptraj(temptraj(:,8)>killsidearm,:)=[];
        end
        
        
        breaks=find(diff(temptraj(:,1))>1); % only use time breaks that last more than a second
        % concatenate the start of each traj, its end, and the origin and
        % destination at that first index
        runepochs=[[temptraj(1); temptraj(breaks+1,1)] [temptraj(breaks,1);...
            temptraj(end,1)] temptraj([1; breaks+1],[4 5])];
        runinds=[[1; breaks+1] [breaks; size(temptraj,1)] temptraj([1; breaks+1],[4 5])];
        %runepochs(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        %runinds(diff(runepochs(:,[1 2]),1,2)<=.5,:)=[]; % remove short epochs
        
        % here would be the place to pick the runs that follow each correct
        % odor sample...
        
        OKruninds=[]; OKrunepochs=[];
        if onlyodorruns
            for rn=1:length(odorStarts)
                % find the next run start
                nextind=find(runepochs(:,1)>odorStarts(rn,1),1,'first');
                OKrunepochs(rn,:)=runepochs(nextind,:);
                OKruninds(rn,:)=runinds(nextind,:);
            end
        else
            OKrunepochs=runepochs;
            OKruninds=runinds;
        end
        
        % identify each run epoch by which type of run it is
        % now get the likelihood of each spike per msec
        runIDs{ile}=ismember(OKruninds(:,[3 4]),trajinds(1,:),'rows');
        runSpkMat=[]; foundRunSpikes={};
        for j=1:length(SuperRat(ses).units)
            [~,runSpkMat{ile}(:,j),~,~,~,foundRunSpikes{:,j}]=event_spikes(SuperRat(ses).units(j).ts,...
                OKrunepochs(:,1),0,OKrunepochs(:,2)-OKrunepochs(:,1));
        end
    end


    thisdecode=mytemplate;
    for i=1:length(region)
        inRegion=cellfun(@(a) contains(a,region{i},'IgnoreCase',true), {SuperRat(ses).units.area});
        if sum(inRegion)>5
            SpkMat=odorSpkMat(:,inRegion & ~badunits)*winsize; % p(spike| 1 ms timebin)
            
            [p_correct,fract_correct_real,fract_correct_null]...
                = nFoldBayesPoisson(SpkMat,odorIDs,5,200);
            
            cv = cvpartition(length(odorIDs), 'kfold',nfolds);
            % xfold crossval
            % this will not work if we ahve cells who are silent during the
            % odor and silent during the run...
            % true decoding
            for k=1:nfolds
                % so dim 1 is going to be test trial, dim 2 is unit and dim 3
                % is right vs left
                odorPriors=mean(SpkMat(odorIDs==1 & cv.training(k),:));
                odorPriors(:,:,2)=mean(SpkMat(odorIDs==0 & cv.training(k),:)); % 1 by nunits by 2 odors
                odorPriors(odorPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
                testmat=SpkMat(cv.test(k),:); % the m trials x n units matrix
                trueOdors=odorIDs(cv.test(k),:); % M trial ids
                
                %p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
                % C  normalizes sum probability to equal 1
                prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
                summat=sum(odorPriors,2); % sum across units (mean_i)
                tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
                probmat=prodmat.*exp(-tau.*summat);
                % sum p across potential odors must ==1 (the c term)
                realProb=probmat./sum(probmat,3);
                [~,decoded]=max(squeeze(realProb),[],2);
                p_success(k)=nanmean(trueOdors==odorNames(decoded)');
            end
            fract_correct_real = mean(p_success);
            
            % do shuffle
            wb=waitbar(0,'Starting Boots');
            for boot = 1:nBoots
                waitbar(boot/nBoots,wb,'Running Boots');
                shuff_odorIDs=odorIDs(randperm(length(odorIDs))); % reorder the matrix so it doesnt match odors
                for k=1:nfolds
                    odorPriors=mean(SpkMat(shuff_odorIDs==1 & cv.training(k),:));
                    odorPriors(:,:,2)=mean(SpkMat(shuff_odorIDs==0 & cv.training(k),:)); % 1 by nunits by 2 odors
                    odorPriors(odorPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
                    testmat=SpkMat(cv.test(k),:); % the m trials x n units matrix
                    trueOdors=shuff_odorIDs(cv.test(k),:); % M trial ids
                    prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
                    summat=sum(odorPriors,2); % sum across units (mean_i)
                    tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
                    probmat=prodmat.*exp(-tau.*summat);
                    % sum p across potential odors must ==1 (the c term)
                    realProb=probmat./sum(probmat,3);
                    [~,decoded_null]=max(squeeze(realProb),[],2);
                    p_success_null(k)=nanmean(trueOdors==odorNames(decoded_null)');
                end
                fract_correct_null(boot) = mean(p_success_null);
            end
            p_correct=1-normcdf(fract_correct_real,nanmean(fract_correct_null),...
                nanstd(fract_correct_null));
            
            
            % now decode runs
            odorPriors=mean(SpkMat(odorIDs==1,:)); % 1,nunits,1
            odorPriors(:,:,2)=mean(SpkMat(odorIDs==0,:)); % 1,nunits,2
            testmat=runSpkMat(:,inRegion & ~badunits); % ntrials x nunits
            trueOdors=runIDs; %
            
            prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
            summat=sum(odorPriors,2);
            prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
            summat=sum(odorPriors,2); % sum across units (mean_i)
            tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
            probmat=prodmat.*exp(-tau.*summat);
            % sum p across potential odors must ==1 (the c term)
            realProb=probmat./sum(probmat,3);
            [~,decoded]=max(squeeze(realProb),[],2);
            fract_correct_run=nanmean(trueOdors==odorNames(decoded)');
            p_run=1-normcdf(fract_correct_run,nanmean(fract_correct_null),...
                nanstd(fract_correct_null));
            
            
            close(wb);
            
            fprintf('%s %d %s trials:%d,  odor decoding correct%.f%%, run decoding correct%.f%%,  null=%.f%% p=%.3f, p%.3f\n',...
                SuperRat(ses).name,SuperRat(ses).daynum,region{i},length(odorIDs),...
                fract_correct_real*100, fract_correct_run*100,...
                nanmean(fract_correct_null)*100,p_correct,p_run);
            
            thisdecode{i+1,2,1}=[fract_correct_real nanmean(fract_correct_null)];
            thisdecode{i+1,3,1}=[fract_correct_run nanmean(fract_correct_null)];
            thisdecode{i+1,2,2}=p_correct; thisdecode{i+1,3,2}=p_run;
            
            
        else

            thisdecode{i+1,2,1}=[nan nan];
            thisdecode{i+1,3,1}=[nan nan];
            thisdecode{i+1,2,2}=nan; thisdecode{i+1,3,2}=nan;
            fprintf('%s %d %s not used, not enough units... \n',...
                SuperRat(ses).name,SuperRat(ses).daynum,region{i});
        end
    end
    SuperRat(ses).BayesDecodePvals=thisdecode;
    % this method of taking the mean run probably doesnt make much sense
    % i think i'd rather do something like a matched timewindow analysis
    
end
warning('on','all');

% so its rows are ca1,pfc, and cols are 1,2,3,4
realProbs{ile,1}=cellfun(@(a) a{2,3,1}(1), {SuperRat.BayesDecodePvals});
nullProbs{ile,1}=cellfun(@(a) a{2,3,1}(2), {SuperRat.BayesDecodePvals});

realProbs{ile,2}=cellfun(@(a) a{3,3,1}(1), {SuperRat.BayesDecodePvals});
nullProbs{ile,2}=cellfun(@(a) a{3,3,1}(2), {SuperRat.BayesDecodePvals});


end
%}