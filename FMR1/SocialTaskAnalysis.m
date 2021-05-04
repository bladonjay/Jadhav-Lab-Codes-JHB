%% generating latin cubes for random animal pairings

% this bit still needs work, I need to make a sequential algorithm that
% checks the followings:
% 1. consec days dont share a first or last pair, or a first or last rat
% too many times
% 2. rats dont have the same partners in the same order more than 2 days in
% a row

%{
ratnames= {'AH1','AH2','AH3','AH6'};
ratairs={};
iterator=1;
for i=1:length(ratnames)-1
    for j=(i+1):length(ratnames)
        ratpairs{iterator}=[ratnames{i} '-' ratnames{j}];
        iterator=iterator+1;
    end
end

sets = [1 2 3 4 5 6;...
        1 5 4 6 3 2;...
        4 1 3 5 2 6;...
        5 3 1 6 2 4;...
        6 5 1 4 3 2]

% build
cubes={}; iterator = 1;
for i=1:size(sets,1)
    thiscube=[];
    for j=1:6
        cubes{i}(j,:)=circshift(1:6,sets(i,j));
        iterator=iterator+1;
    end
    
end

% now from each day... this needs serious work...
iterator=1;
randompairs={};
for i=1:6 % six perms for weeks
    % pull one random col from each cube successively until you've filled
    % for good
    for j=1:size(sets,1) % each day
        for k=1:6 % each pairing in a day
            randompairs{i}{j,k}=ratpairs{cubes{j}(k,i)};
        end
    end
end
%}

%% aggregating all the pairings...

try
    load SocialSessions.m
catch
    
    opts = delimitedTextImportOptions("NumVariables", 2);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = "#";
    
    % Specify column names and types
    opts.VariableNames = ["VarName1", "Var2"];
    opts.SelectedVariableNames = "VarName1";
    opts.VariableTypes = ["string", "string"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, ["VarName1", "Var2"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["VarName1", "Var2"], "EmptyFieldRule", "auto");
    
    mydir=uigetdir;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% go through and add metadata %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GroupPrefix='20';
    idlength=2;
    
    
    filelist = dir(fullfile(mydir, '**\*.*'));  %get list of files and folders in any subfolder
    %filelist=getAllFiles(mydir,'.stateScriptLog');
    okfiles=cellfun(@(a) contains(a,'stateScriptLog'),{filelist.name});
    filelist = filelist(okfiles);
    
    % now get the rat and date
    for i=1:length(filelist)
        try
            sessIDX=strfind(filelist(i).name,'(');
            datestr=filelist(i).name(4:sessIDX-1);
            namestr=filelist(i).name(sessIDX:strfind(filelist(i).name,'.')-1);
            nameIDX=strfind(filelist(i).name(sessIDX:end),'-');
            numIDx=isstrprop(namestr,'digit');
            
            %filelist(i).rundate=filelist(i).folder(find(filelist(i).folder=='\',1,'last')+1:end);
            filelist(i).rundate=datestr;
            filelist(i).datenum=datenum(datestr);
            
            filelist(i).runnum=namestr(find(numIDx,1,'first'));
            
            lastdigits=find(diff(numIDx)==-1);
            
            filelist(i).ratnum(1)=str2double(namestr(lastdigits(2)));
            filelist(i).ratnum(2)=str2double(namestr(lastdigits(3)));
            filelist(i).ratCohort=namestr(lastdigits(2)-2 : lastdigits(2)-1);
        end

    end
    
    % now sort, first by session, then by date then by rat, this will make rat
    % the top sorting category
    
    % want it as a struct or a table?
    filelist=filelist(cellfun(@(a) ~isempty(a), {filelist.ratnum}));
    % we can interchange tables and structs easily
    %ratinfo=sortrows(struct2table(filelist),{'datenum','sessnum'});
    
    % now for ease here im going to put it into a struct
    %ratinfo=table2struct(ratinfo);
    
    
end
%%


% okay now we can interpret these data...
% the question is what do we want to know?

% the essential question is whether animals are actually following
% eachother

% I can think of a few ways to examine this
% 1. how many arm arrivals result in a reward?
% 2. do the departure times after reward correlate between animals?
% 3. who gets to the rewarded arms first? is this consistent?
% 4. are their arrival times correlated? (wont be useful with new code)
% 5. how many rewards per minute?
% 6. is a rat more likely to enter an arm his friend is in vs an empty one?
        % this will cahnge based on who leads or lags

% control questions
% 1. do animals prefer certain buddies (do their visits/transitions change
% with buddy?
% 2. does performance change across the day?

%%



%%



verbose=0;
ratinfo=filelist;
ratinfo=ratinfo(cellfun(@(a) ~contains(lower(a),'z'),{ratinfo.name}));
%%
for i=1:length(ratinfo)
    fprintf(' \n \n');
    fprintf('Running cohort %s sess %d  with %s%d and %s%d \n',ratinfo(i).ratCohort,i,...
        ratinfo(i).ratCohort, ratinfo(i).ratnum(1), ratinfo(i).ratCohort, ratinfo(i).ratnum(2));
    DataFile = readtable(fullfile(ratinfo(i).folder,ratinfo(i).name), opts);
    % Convert to output type
    DataRaw = table2cell(DataFile);
    numIdx = cellfun(@(x) ~isempty(x{1}), DataRaw(:,1));
    DataTips = DataRaw(numIdx,1);
    
    % convert to char and split out
    DataTips2=cellfun(@(a) char(a{1}), DataTips, 'UniformOutput',false);
    DataAll=cellfun(@(a) a(1:find(a==' ',1,'first')), DataTips2,'UniformOutput',false);
    DataAll(:,2)=cellfun(@(a) a(find(a==' ',1,'first')+1:end), DataTips2,'UniformOutput',false);
    
    ledger=DataAll;
    myevents={};
    % here we must swap for parseSocialEvents I think...
    % we split up into two sets of events, and we debounce repeated hits on
    % the same well (if the last end and the next start are within a
    % second, its the same sample
    debounce=1; % 1 second debounce...
    [myevents{1},myevents{2}] = parseSocialEvents(ledger,debounce);

    if verbose
        % make a huge plot here
        % on left y maybe a moving average of total success rate
        figure;
        for rt=1:length(myevents)
        for it=1:length(myevents{rt})
            if myevents{rt}(it,5)==0
                plot([myevents{rt}(it,1) myevents{rt}(it,2)],repmat(myevents{rt}(it,3),1,2)+3*rt,'k-','LineWidth',3);
                hold on
            else
                plot([myevents{rt}(it,1) myevents{rt}(it,2)],repmat(myevents{rt}(it,3),1,2)+3*rt,'r-','LineWidth',3);
                hold on
            end
        end
        end
    end
    
    % tables are easier to read...
    for rt=1:2
        ratinfo(i).ratsamples{rt}=table(myevents{rt}(:,1),myevents{rt}(:,2),...
            myevents{rt}(:,3),myevents{rt}(:,4),myevents{rt}(:,5),...
            myevents{rt}(:,6), myevents{rt}(:,7),'VariableNames',...
            {'start','end','thiswell','lastwell','Reward','match','Goal Well'});
    
        % now report
        ctrls=myevents{rt}(diff(myevents{rt}(:,[3 4]),1,2)~=0,:);
        fprintf('Today rat %d, %.2f%% arm transitions were to an occupied arm, and %.2f%% were rewarded \n',...
            ratinfo(i).ratnum(rt), nanmean(ctrls(:,6))*100,nanmean(ctrls(:,5))*100); % basically when the diff==0
    
    end
    ratinfo(i).datenum=datenum(ratinfo(i).date);

end


%% okay lets really analyze some of these questions

% first, i want to see whether the animals actually tend to go to an
% occupied arm more often than a random one
    
% to do this, we will need to restrict analyses to times when the other rat
% is at one of the other two arms, and then our rat has a choice between a
% and b.

% need to get for each arm transition, the most recent arm the other rat
% occupied, and when he got there and when/if he left.


% so lets try one rat here
% for each session, gather if this rat is the rat, if it is, then its id
% goes in the first column, and its buddy in the peer column
% then for each of his transitions, we can fin

rattable=struct2table(ratinfo);
sessRats=cat(1,rattable.ratnum');
allRatNames=unique(sessRats(:));
RatAll=struct('names',allRatNames);
%ratTable=table(allRatNames,'VariableNames',{'Rat Name'});
myCmap=lines(4);

doStateSpace=0;

for ra=1:length(allRatNames)
    iterator=1;
    for ses=1:length(ratinfo)
        
        ratmatch=allRatNames(ra)==sessRats(:,ses);
        if any(ratmatch)
            myevents=ratinfo(ses).ratsamples{ratmatch};
            hisevents=ratinfo(ses).ratsamples{~ratmatch};
            
            % zero out everything to the start!
            firstevent=min([myevents.start; hisevents.start]);
            myevents.start=myevents.start-firstevent;
            myevents.end=myevents.end-firstevent;
            hisevents.start=hisevents.start-firstevent;
            hisevents.end=hisevents.end-firstevent;

            
            
            
            if iterator==1, myRaw=myevents; hisRaw=hisevents;
            else
                if myevents.start(1)<myRaw.end(end)
                    lastevent=max([myRaw.end; hisRaw.end])+rand*10;
                    myevents.start=myevents.start+lastevent;
                    myevents.end=myevents.end+lastevent;
                    hisevents.start=hisevents.start+lastevent;
                    hisevents.end=hisevents.end+lastevent;
                end
                
                myRaw=[myRaw; myevents]; hisRaw=[hisRaw; hisevents];
            end
            
            
            myevents.hiswell=zeros(height(myevents),1);
            myevents.arrival=zeros(height(myevents),1);
            myevents.departure=zeros(height(myevents),1);
            myevents.sess=zeros(height(myevents),1);
            for i=1:height(myevents)
                % for each event, find the most recent event of his
                hiseventID=find(hisevents.start<myevents.start(i),1,'last');
                if ~isempty(hiseventID)
                    myevents.hiswell(i)=hisevents.thiswell(hiseventID);
                    myevents.arrival(i)=hisevents.start(hiseventID);
                    myevents.departure(i)=hisevents.end(hiseventID);
                end
            end
            myevents.sess(:)=iterator;
            if iterator==1 
                allMyEvents=myevents;
            else
                allMyEvents=[allMyEvents; myevents];
            end
            iterator=iterator+1;
        end
    end
    
    % okay first scrub all trials where hiswell is my lastwell
    allMyEvents(allMyEvents.thiswell==allMyEvents.lastwell,:)=[]; % kill return visits
    allMyEvents(allMyEvents.hiswell==allMyEvents.lastwell,:)=[]; % kill  when were at the same well
    
    % take all the arm transitions now across time
    RatAll(ra).data=allMyEvents;
    RatAll(ra).myRaw=myRaw;
    RatAll(ra).hisRaw=hisRaw;
    
    allMyEvents(allMyEvents.departure<allMyEvents.start-2,:)=[]; % kill when he left that well too long ago
    % lets do a cumsum learning curve
    matches=allMyEvents.thiswell==allMyEvents.hiswell;
    if doStateSpace
        [bmode,b05,b95,pmatrix,wintr] = CalcStateSpacePerformance(matches', .5,0);
        %plot(cumsum(matches)./(1:height(allMyEvents))');
        pmatrix(1:10)=nan;
        wintr=find(pmatrix<.05,1,'first');
        if ~isnan(wintr)
            pl(5)=plot(1:wintr,bmode(1:wintr),':','color',myCmap(ra,:),'LineWidth',2);
            hold on;
            pl(ra)=plot(wintr:length(bmode),bmode(wintr:end),'-','color',myCmap(ra,:),'LineWidth',2);
        else
            pl(ra)=plot([0 0],[1 1],'-','color',myCmap(ra,:),'LineWidth',2);
            plot(bmode,':','color',myCmap(ra,:),'LineWidth',2);
        end
        
        
        RatAll(ra).bmode=bmode;
        RatAll(ra).pmatrix=pmatrix;
    end
    % now kill all where he leaves more than 2 seconds before we get there
end

pl(6)=plot([1 800],[.5 .5],'k');
xlim([0 650]);
legend(pl,{'Rat 1','Rat 2','Rat 3','Rat 4','Before Social Preference','Baseline'});
ylabel(sprintf('Likelihood of choosing a \n peer-occupied arm'));
xlabel('Reward Visit Number');
% so result here is that 2/4 animals learned to go to arms where their
% friend was occupying in order to get rewards faster.

% the basis is that when given the choice, 2/4 of these animals chose the
% arm their friend was at over the arm they were at

%% next question...

% howabout their efficiency, does it get better?  what is baseline
% efficiency?  my guess is something like likelihood of other guy being at
% an arm, * 1/3 for three arms.. so what %% of time is partner at an arm,
% multiply that by 1/3 and you get baseline chances of finding him...


%% how about whether they move together?

% basically do they alternate wells in close proximity?
% the question i guess is what is the distribution of latencies between
% animal arrival times, and is this greater than chance?

% I think the way to do this is to take the xcorr of these arrival times at
% say 1-30 second lags, and then see if there are any peaks... then run
% this same analysis after circshifting one.

% the other question is whether this actually swaps from time to time...
% maybe one animal leads and so theres a significant bump behind him, and
% then sometimes he follows, so the bump is ahead of him.  It may be a good
% idea to do this within each session.  The best result would be that the
% distribution is multimodal, and each animal has their favorite lag



% one thing to do would be to remove all arrivals immediately following a
% reward, because then they;re locked to whenever the rats finish eating
figure;
for ra=1:4
    
    binsize=1/5; maxlag=120;
    timeshifts=(-maxlag:maxlag)*binsize;
    mytrans=RatAll(ra).myRaw(RatAll(ra).myRaw.thiswell~=RatAll(ra).myRaw.lastwell,:);
    myarrivals=accumarray(round(mytrans.start./binsize)+1,1);
    
    histrans=RatAll(ra).hisRaw(RatAll(ra).hisRaw.thiswell~=RatAll(ra).hisRaw.lastwell,:);
    hisarrivals=accumarray(round(histrans.start./binsize)+1,1);
    
    
    if length(hisarrivals)<length(myarrivals)
        hisarrivals((end+1):length(myarrivals))=0;
    elseif length(hisarrivals)>length(myarrivals)
        myarrivals(end:length(hisarrivals))=0;
    end
    myarrivals=SmoothMat2(myarrivals,[10 10],2);
    hisarrivals=SmoothMat2(hisarrivals,[10 10],2);
    
    % right direction means shift second input foreward
    temp=xcorr(myarrivals,hisarrivals,maxlag);
    subplot(2,2,ra);
    b=bar(timeshifts,temp,1,'EdgeColor','none','FaceColor',myCmap(ra,:));
    bootcorr=[];
    % now for a bootstrap
    for i=1:200
        myboot=circshift(myarrivals,randi(length(myarrivals)));
        bootcorr(i,:)=xcorr(myboot,hisarrivals,maxlag);
    end
    
    hold on;
    b(2)=plot(timeshifts,prctile(bootcorr,5),'r');
    plot(timeshifts,prctile(bootcorr,95),'r');
    ylim([min(prctile(bootcorr,5))*.8,max(prctile(bootcorr,95))*1.2]);
    plot([0 0],[0 max(temp)*1.2],'k');
    title(sprintf('Rat %d',ra));
    switch ra
        case 1
            ylabel('Correlation (au)');
        case 3
            xlabel(sprintf(' Leads Partner <--> Follows partner \n Temporal Offset in seconds'));
            ylabel('Correlation (au)');
        case 4
            xlabel(sprintf(' Leads Partner <--> Follows partner \n Temporal Offset in seconds'));
    end
end
legend(b,'Real Data','95% CI');

figure;
for ra=1:4
    % remove rewards
    noReward=RatAll(ra).data(2:end,:);
    killtr=find(noReward.Reward==1)+1;
    killtr=killtr(killtr<=(height(noReward)-2));
    noReward(killtr,:)=[];
    
    % remove any zeros
    killtr=noReward.start==0 | noReward.arrival==0;
    noReward(killtr,:)=[];
    
    % accumarray
    hisarrivals=accumarray(round(noReward.arrival),1);
    myarrivals=accumarray(round(noReward.start),1);
    
    % match array lengths
    if length(hisarrivals)<length(myarrivals)
        hisarrivals((end+1):length(myarrivals))=0;
    elseif length(hisarrivals)>length(myarrivals)
        myarrivals(end:length(hisarrivals))=0;
    end
    maxlag=120; % 2 minuites
    
    temp=xcorr(myarrivals,hisarrivals,maxlag);
    subplot(2,2,ra);
    bar(-maxlag:maxlag,temp,'FaceColor',myCmap(ra,:));
    bootcorr=[];
    % now for a bootstrap
    for i=1:1000
        myboot=circshift(myarrivals,randi(length(hisarrivals)));
        bootcorr(i,:)=xcorr(myboot,hisarrivals,maxlag);
    end
    
    hold on;
    plot(-maxlag:maxlag,prctile(bootcorr,1),'r');
    plot(-maxlag:maxlag,prctile(bootcorr,99),'r');
    plot([0 0],[0 max(temp)*1.2],'k');
    title(RatAll(ra).names);
    switch ra
        case 1
            ylabel('Correlation (au)');
        case 3
            xlabel('Temporal Offset in seconds');
            ylabel('Correlation (au)');
        case 4
            xlabel(sprintf(' Leads Partner <--> Follows partner \n Temporal Offset in seconds'));
    end
    axis tight
end
sgtitle('Each rat alternates arms close in time with his partner');














%%
% howabout a number of wins per arm visits for each rat, and for each pair
% type


% looks like the animals strategy oscillates with who their partner is
figure;
for i=1:4
    plot(SmoothMat2(RatAll(i).myRaw.match,[1 250],100),'--','Color',myCmap(i,:));
    hold on;
    plot(SmoothMat2(RatAll(i).myRaw.Reward(:),[1 250],100),'-','Color',myCmap(i,:));
end
    

%% howabout the efficiency of each session?

% so i guess the question is this, per arm visit, are the control pairs
% getting more rewards?


% THIS IS QUICK AND DIRTY, AND ALSO CONSIDERS THE NUMBER OF ARM TRANSITIONS
% A RAT DOES EVEN WHEN HIS PARTNER IS AT THE WELL HE IS DEPARTING


% first quantify number of arm visits per animal
ctrls=[]; fxpair=[]; combos=[];
ctrlsraw={}; fxpairraw={}; comboraw={};

[sessnums,ia,ic]=unique(cellfun(@(a) datenum(a), {ratinfo.rundate}));
for i=1:length(sessnums)
    daysess=ratinfo(ic==i);
    cumct=1;
    for k=1:length(daysess)       
    daysess(k).ratsamples{1}.Properties.VariableNames{6} = 'match';
    daysess(k).ratsamples{2}.Properties.VariableNames{6} = 'match';
    if any(daysess(k).ratnum==1) &&...
            any(daysess(k).ratnum==4) &&...
            str2double(daysess(k).runnum)<7
        % quantify # arm transitions for each rat
        ctrls(i,1)=sum(diff(daysess(k).ratsamples{1}.thiswell)~=0);
        ctrls(i,2)=sum(diff(daysess(k).ratsamples{2}.thiswell)~=0);
        % and now the number of matched samples
        % its the same for each rat because both rats sample when its a
        % match
        ctrls(i,3)=sum(daysess(k).ratsamples{2}.match);
        ctrlsraw{i,1}=daysess(k).ratsamples{1}.match(diff(daysess(k).ratsamples{1}.thiswell)~=0);
        ctrlsraw{i,2}=daysess(k).ratsamples{2}.match(diff(daysess(k).ratsamples{2}.thiswell)~=0);
    elseif any(daysess(k).ratnum==2) &&...
            any(daysess(k).ratnum==3) && ...
            str2double(daysess(k).runnum)<7
        fxpair(i,1)=sum(diff(daysess(k).ratsamples{1}.thiswell)~=0);
        fxpair(i,2)=sum(diff(daysess(k).ratsamples{2}.thiswell)~=0);
        % and now the number of rewards
        fxpair(i,3)=sum(daysess(k).ratsamples{2}.match);
        fxpairraw{i,1}=daysess(k).ratsamples{1}.match(diff(daysess(k).ratsamples{1}.thiswell)~=0);
        fxpairraw{i,2}=daysess(k).ratsamples{2}.match(diff(daysess(k).ratsamples{2}.thiswell)~=0);

    elseif str2double(daysess(k).runnum)<7
        combos(i,1,cumct)=sum(diff(daysess(k).ratsamples{1}.thiswell)~=0);
        combos(i,2,cumct)=sum(diff(daysess(k).ratsamples{2}.thiswell)~=0);
        % and now the number of rewards
        combos(i,3,cumct)=sum(daysess(k).ratsamples{2}.match);
        fxpairraw{i,1,cumct}=daysess(k).ratsamples{1}.match(diff(daysess(k).ratsamples{1}.thiswell)~=0);
        fxpairraw{i,2,cumct}=daysess(k).ratsamples{2}.match(diff(daysess(k).ratsamples{2}.thiswell)~=0);
        cumct=cumct+1;
    end
    end
end

ctrls([12 22],:)=[];
fxpair([12 22],:)=[];
combos([12 22],:,:)=[];
combomeans=nanmean(combos,3);
figure;
subplot(1,2,1);
% %% arm transitions go to friends arm
plot(ctrls(:,3)./mean(ctrls(:,1:2),2)); hold on;
plot(fxpair(:,3)./mean(fxpair(:,1:2),2));
plot(combomeans(:,3)./mean(combomeans(:,1:2),2));

title(sprintf('ranksum test p= %.2e',...
    ranksum(ctrls(:,3)./max(ctrls(:,1:2),[],2),fxpair(:,3)./max(fxpair(:,1:2),[],2))));
xlabel('Session number');
ylabel(sprintf('Proportion of arm transitions \n to peer-occupied arm'));
legend('WT-WT pair','FX-FX pair','FX-WT pair');
% i wonder if their armtransitions per minute were higher.
ctrls=[]; fxpair=[];
cumct=[1 1];
for i=1:length(ratinfo)
    if any(ratinfo(i).ratnum==1) &&...
            any(ratinfo(i).ratnum==4) &&...
            str2double(ratinfo(i).runnum)<7
        % session start, session end
        ctrls(cumct(1),1)=ratinfo(i).ratsamples{1}.start(1);
        ctrls(cumct(1),2)=ratinfo(i).ratsamples{1}.end(end);
        % and now the number of transitions that day
        ctrls(cumct(1),3)=sum(diff(ratinfo(i).ratsamples{1}.thiswell)~=0)+...
            sum(diff(ratinfo(i).ratsamples{2}.thiswell)~=0);
        cumct(1)=cumct(1)+1;
    elseif any(ratinfo(i).ratnum==2) &&...
            any(ratinfo(i).ratnum==3) && ...
            str2double(ratinfo(i).runnum)<7
        fxpair(cumct(2),1)=ratinfo(i).ratsamples{1}.start(1);
        fxpair(cumct(2),2)=ratinfo(i).ratsamples{1}.end(end);
        % and now the number of rewards
        fxpair(cumct(2),3)=sum(diff(ratinfo(i).ratsamples{1}.thiswell)~=0)+...
            sum(diff(ratinfo(i).ratsamples{2}.thiswell)~=0);
        cumct(2)=cumct(2)+1;
    end
end

ctrls([12 17 23],:)=[];
fxpair([12 17 23],:)=[];
subplot(1,2,2);
plot(ctrls(:,3)./(ctrls(:,2)-ctrls(:,1)).*60); hold on;
plot(fxpair(:,3)./(fxpair(:,2)-fxpair(:,1)).*60);
xlabel('Session number');
ylabel('Number of arm transitions per minute');
legend('WT-WT pair','FX-FX pair');

title(sprintf('ranksum test p= %.2e',...
    ranksum(ctrls(:,3)./max(ctrls(:,1:2),[],2),fxpair(:,3)./max(fxpair(:,1:2),[],2))));
%%

for i=1:length(ratinfo)
    ratinfo(i).sessnum=str2double(ratinfo(i).runnum);
    ratinfo(i).datenum=datenum(ratinfo(i).date);
end