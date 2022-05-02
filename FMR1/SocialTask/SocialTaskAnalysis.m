 %% generating latin cubes for random animal pairings

% this bit still needs work, I need to make a sequential algorithm that
% checks the followings:
% 1. consec days dont share a first or last pair, or a first or last rat
% too many times
% 2. rats dont have the same partners in the same order more than 2 days in
% a row



%
%
% for spike data use gatherMS_to_FFunits in parent folder
%
%

edit gatherMS_to_FFunits.m

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
% first get your excel doc
load('G:\SocialData\Behavior only\SocialCohorts1-3.mat')

try
    %load SocialCohorts1-3.m
    load('G:\SocialData\Behavior only\SocialCohorts1-3.mat')

catch
    
    mydir=uigetdir;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% go through and add metadata %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 201-202,  XFZ 1-2, XFB1-XFB2
    GroupPrefix={'20','XFB','XFZ','XFE','XFG'};

    
    
    fileList = dir(fullfile(mydir, '**\*.*'));  %get list of files and folders in any subfolder

    okfiles=cellfun(@(a) contains(a,'stateScriptLog'),{fileList.name});
    fileList = fileList(okfiles);
    
    okfiles=ones(length(fileList),1);
    
    % now get the rat and date
    for i=1:length(fileList)
        try
            startIDX=strfind(fileList(i).name,'(');
            endIDX=strfind(fileList(i).name,')');
            
            % kill test files
            if contains(fileList(i).name(startIDX:endIDX),'test')
                okfiles(i)=0;
            else
                
                % pull dates
                %datestr=fileList(i).name(4:startIDX-1);
                realdate=fileList(i).date(1:find(fileList(i).date==' ')-1);
                fileList(i).runtime=fileList(i).date(find(fileList(i).date==' ')+1:end);
                fileList(i).datenum=datenum(realdate);
                
                namestr=fileList(i).name(startIDX+1:endIDX-1);
                namestart=find(~isstrprop(namestr,'digit'),1,'first');
                fileList(i).runnum=namestr(1:namestart-1);
                if strcmpi(namestr(namestart+1:namestart+2),'20')
                    fileList(i).cohortnum=1;
                    fileList(i).cohortname='20';
                    fileList(i).ratnames={namestr(namestart+3),namestr(namestart+7)};
                    fileList(i).ratnums=[str2num(namestr(namestart+3)) str2num(namestr(namestart+7))];
                elseif strcmpi(namestr(namestart+1:namestart+3),'XFB')
                    fileList(i).cohortnum=2;
                    fileList(i).cohortname='XFB';
                    fileList(i).ratnames={namestr(namestart+4), namestr(namestart+9)};
                    fileList(i).ratnums=[str2num(namestr(namestart+4)) str2num(namestr(namestart+9))];
                elseif strcmpi(namestr(namestart+1:namestart+3),'XFZ')
                    %6 and 8?
                    fileList(i).cohortnum=4;
                    fileList(i).cohortname='XFZ';
                    fileList(i).ratnames={namestr(namestart+6), namestr(namestart+8)};
                    fileList(i).ratnums=[str2num(namestr(namestart+6)) str2num(namestr(namestart+8))];
                elseif strcmpi(namestr(namestart+1:namestart+3),'XFE')
                    %6 and 8?
                    fileList(i).cohortnum=5;
                    fileList(i).cohortname='XFE';
                    fileList(i).ratnames={namestr(namestart+4), namestr(namestart+9)};
                    fileList(i).ratnums=[str2num(namestr(namestart+4)) str2num(namestr(namestart+9))];
                elseif strcmpi(namestr(namestart+1:namestart+3),'XFG')
                    %6 and 8?
                    fileList(i).cohortnum=3;
                    fileList(i).cohortname='XFG';
                    fileList(i).ratnames={namestr(namestart+4), namestr(namestart+9)};
                    fileList(i).ratnums=[str2num(namestr(namestart+4)) str2num(namestr(namestart+9))];
                else
                    fprintf('Cant parse this, the filename is... \n     %s \n',...
                        fileList(i).name);
                    okfiles(i)=0;
                end
            end
        catch
            okfiles(i)=0;
        end
        
    end
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
% set options for data import and only import data we want
if ~isfield(ratinfo,'ratsamples')
    fileList=fileList(okfiles==1);
    %fileList=fileList(sum([fileList.cohortnum]==[1; 3])'>0); % right now its all but you'll want to pull only certain groups
    opts=setStateScriptOpts();
    
    verbose=0;
    ratinfo=fileList;
    ratinfo=ratinfo(cellfun(@(a) ~contains(lower(a),'z'),{ratinfo.name}));
    ratinfo=ratinfo(cellfun(@(a) a>0, {ratinfo.bytes}));
    %
    for i=1:length(ratinfo)
        fprintf(' \n \n');
        fprintf('Running cohort %d sess %d  with %s%d and %s%d \n',ratinfo(i).cohortnum,i,...
            ratinfo(i).cohortname, ratinfo(i).ratnums(1), ratinfo(i).cohortname, ratinfo(i).ratnums(2));
        DataFile = readtable(fullfile(ratinfo(i).folder,ratinfo(i).name), opts);
        % Convert to output type
        
        myevents={};
        % here we must swap for parseSocialEvents I think...
        % we split up into two sets of events, and we debounce repeated hits on
        % the same well (if the last end and the next start are within a
        % second, its the same sample
        debounce=1; % 1 second debounce...
        [myevents{1},myevents{2}] = parseSocialEvents(DataFile,debounce);
        
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
        fprintf('Tot rewards for pair %s%d & %s%d was %d\n',ratinfo(i).cohortname,...
            ratinfo(i).ratnums(1), ratinfo(i).cohortname, ratinfo(i).ratnums(2), sum(myevents{1}(:,5)));
        for rt=1:2
            ratinfo(i).ratsamples{rt}=table(myevents{rt}(:,1),myevents{rt}(:,2),...
                myevents{rt}(:,3),myevents{rt}(:,4),myevents{rt}(:,5),...
                myevents{rt}(:,6), myevents{rt}(:,7),'VariableNames',...
                {'start','end','thiswell','lastwell','Reward','match','Goal Well'});
            
            % now report
            firsthit=myevents{rt}(diff(myevents{rt}(:,[3 4]),1,2)~=0,:);
            fprintf('     rat %d, %.2f%% transitions were to occupied arm, %.2f%% were rewarded \n',...
                ratinfo(i).ratnums(rt), nanmean(firsthit(:,6))*100,nanmean(firsthit(:,5))*100); % basically when the diff==0
            
        end
        
    end
end

%% okay lets really analyze some of these questions
%{
first, i want to see whether the animals actually tend to go to an
occupied arm more often than a random one
    
to do this, we will need to restrict analyses to times when the other rat
is at one of the other two arms, and then our rat has a choice between a
and b.

need to get for each arm transition, the most recent arm the other rat
occupied, and when he got there and when/if he left.


so lets try one rat here
for each session, gather if this rat is the rat, if it is, then its id
goes in the first column, and its buddy in the peer column
then for each of his transitions, we can fin
%}
%% statespace performance by rat (it oscillates and is really dirty

% this needs a lot of work
%{

%ratTable=table(allRatNames,'VariableNames',{'Rat Name'});
myCmap=lines(4);

doStateSpace=0;

mypairs=[1 2; 4 3];
for cohort=1
    thisCohort=find((cell2mat({ratinfo.cohortnum})==cohort));
    
for ra=1:2
    itr=1;
    for ses=1:length(thisCohort)
        
        %ratmatch=mypairs(:,ra)==ratinfo(thisCohort(ses)).ratnum;
        %if sum(ratmatch(:))==2
            myevents=ratinfo(thisCohort(ses)).ratsamples{1};
            hisevents=ratinfo(thisCohort(ses)).ratsamples{2};
            
            % zero out everything to the start!
            firstevent=min([myevents.start; hisevents.start]);
            myevents.start=myevents.start-firstevent;
            myevents.end=myevents.end-firstevent;
            hisevents.start=hisevents.start-firstevent;
            hisevents.end=hisevents.end-firstevent;


            if itr==1, myRaw=myevents; hisRaw=hisevents;
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
            myevents.sess(:)=itr;
            if itr==1 
                allMyEvents=myevents;
            else
                allMyEvents=[allMyEvents; myevents];
            end
            itr=itr+1;
        %end
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
end
% so result here is that 2/4 animals learned to go to arms where their
% friend was occupying in order to get rewards faster.

% the basis is that when given the choice, 2/4 of these animals chose the
% arm their friend was at over the arm they were at
%% statespace performance by pair of rats

%{
ratpairs=[1 4; 2 3];
itr=[1 1];
for i=length(ratinfo)
    % first control pair
    if any(ratinfo(i).ratnum==ratpairs(1)) && any(ratinfo(i).ratnum==ratpairs(1,2))
        
        allsamples{1}=ratinfo(i).ratsamples{1};
        allsamples{2}=ratinfo(i).ratsamples{2};
        
        itr=2;
        for k=1:2
            while itr<=height(allsamples{k})
                if allsamples{k}.start(itr)<allsamples{k}.end(itr-1)+5
                    allsamples{k}.end(itr-1)=allsamples{k}.end(itr);
                    allsamples{k}.Reward(itr-1)=allsamples{k}.Reward(itr) || allsamples{k}.Reward(itr-1);
                    allsamples{k}.match(itr-1)=allsamples{k}.match(itr) || allsamples{k}.match(itr-1);
                    allsamples{k}(itr,:)=[];
                else
                    itr=itr+1;
                end
            end
        end
        allevents=matchSocialSamples(allsamples{1},allsamples{2});
        % now look at when the arms transition
        entries=allevents(allevents.in1out0==1,:);
        % for each entry we want to know if its a new well entry
        % need to comb over

%}

%}
%% next question...

% howabout their efficiency, does it get better?  what is baseline
% efficiency?  my guess is something like likelihood of other guy being at
% an arm, * 1/3 for three arms.. so what %% of time is partner at an arm,
% multiply that by 1/3 and you get baseline chances of finding him...


edit ratBehaviorSynchrony


%%
% for each session, keep track of your peers well at all times:

for i=1:length(ratinfo)
    for tr=1:height(ratinfo(i).ratsamples{1}) % for each trial
        % find his last well poke
        hislastpoke=find(ratinfo(i).ratsamples{2}.start<ratinfo(i).ratsamples{1}.start(tr),1,'last');
        if ~isempty(hislastpoke)
            ratinfo(i).ratsamples{1}.hiswell(tr)=ratinfo(i).ratsamples{2}.thiswell(hislastpoke);
        else
            ratinfo(i).ratsamples{1}.hiswell(tr)=nan;
        end
    end
    
    for tr=1:height(ratinfo(i).ratsamples{2})
        hislastpoke=find(ratinfo(i).ratsamples{1}.start<ratinfo(i).ratsamples{2}.start(tr),1,'last');
        if ~isempty(hislastpoke)
            ratinfo(i).ratsamples{2}.hiswell(tr)=ratinfo(i).ratsamples{1}.thiswell(hislastpoke);
        else
            ratinfo(i).ratsamples{2}.hiswell(tr)=nan;
        end
    end
    
end

%% howabout the efficiency of each session?


runtots=0;
% so i guess the question is this, per arm visit, are the control pairs
% best?
% cohort 1 (20x) 1 and 4 are controls
% cohort 2 (XFB) 1 and 2 are controls
% cohort 3 (XFG) 1 and 2 are controls

genotypetable=table([1 4; 1 2; 3 4],[2 3; 3 4; 1 2], 'VariableNames',{'controls','fx'});
figure;
allcohorts=[]; rundates=[];
iterator=1;
for cohort=[1 2 3]
    cohortinfo=ratinfo(cell2mat({ratinfo.cohortnum})==cohort);
    
    % first quantify number of arm visits per animal
    firsthit=[]; fxpair=[]; combos=[];
    
    varnames={'perfA','perfB','totA','totB','totMatches','totRewards'};
    ctrlsraw=table([],[],[],[],[],[],'VariableNames',varnames);
    fxpairraw=ctrlsraw;
    combosraw=table([],[],[],[],[],[],'VariableNames',varnames);
    
    % get unique days, so you can get a day average
    [sessnums,ia,ic]=unique(cellfun(@(a) datenum(a), {cohortinfo.datenum}));
    %rundates=[rundates; [datestr(sessnums) ones(length(sessnums),1)*cohort]];
    for i=1:length(sessnums)
        daysess=cohortinfo(ic==i);
        cumct=1;
        % for each day you want an average of wt wt, wt fx, and fx fx for
        % wins, candidates and the ratio of the two..
        % end matrix will be a m day by 5 by n sess/day matrix
        for k=1:length(daysess)
            
            % controls (no fx?) only one sess
            if ~any(intersect(genotypetable.fx(cohort,:),daysess(k).ratnums)) && ...
                    str2double(daysess(k).runnum)<7
                % quantify # arm transitions for each rat
                for rt=1:2
                    % my last well ~= my current well and my last well ~=
                    % his current well
                    candidates=daysess(k).ratsamples{rt}.thiswell~=daysess(k).ratsamples{rt}.lastwell & ...
                        daysess(k).ratsamples{rt}.hiswell~=daysess(k).ratsamples{rt}.lastwell;
                    % my last well ~= my current well and my current well
                    % == his current well
                    wins=daysess(k).ratsamples{rt}.thiswell~=daysess(k).ratsamples{rt}.lastwell & ...
                        daysess(k).ratsamples{rt}.hiswell==daysess(k).ratsamples{rt}.thiswell;
                    % % of candidates are wins
                    ctrlsraw{i,rt}=sum(wins)/sum(candidates); % get performance rate
                    ctrlsraw{i,rt+2}=sum(candidates); % get total moves
                    %ctrlsraw.rat1tries(i)=sum(candidates);
                    %ctrlsraw.rat1wins(i)=sum(wins);
                end
                % rews / tot trans
                ctrlsraw{i,5}=sum(daysess(k).ratsamples{rt}.match)/sum(ctrlsraw{i,3:4}); % get matches per total moves
                ctrlsraw{i,6}=sum(daysess(k).ratsamples{rt}.Reward)/sum(ctrlsraw{i,3:4}); % get Rewards per total moves
                
                
            % fx only
            elseif ~any(intersect(genotypetable.controls(cohort,:),daysess(k).ratnums)) && ...
                    str2double(daysess(k).runnum)<7
                % quantify # arm transitions for each rat
                for rt=1:2
                    % my last well ~= my current well and my last well ~=
                    % his current well
                    candidates=daysess(k).ratsamples{rt}.thiswell~=daysess(k).ratsamples{rt}.lastwell & ...
                        daysess(k).ratsamples{rt}.hiswell~=daysess(k).ratsamples{rt}.lastwell;
                    % my last well ~= my current well and my current well
                    % == his current well
                    wins=daysess(k).ratsamples{rt}.thiswell~=daysess(k).ratsamples{rt}.lastwell & ...
                        daysess(k).ratsamples{rt}.hiswell==daysess(k).ratsamples{rt}.thiswell;
                    % % of candidates are wins
                    fxpairraw{i,rt}=sum(wins)/sum(candidates); % get performance rate
                    fxpairraw{i,rt+2}=sum(candidates); % get total moves
                end
                % rews / tot trans
                fxpairraw{i,5}=sum(daysess(k).ratsamples{rt}.match)/sum(fxpairraw{i,3:4}); % get matches per total moves
                fxpairraw{i,6}=sum(daysess(k).ratsamples{rt}.Reward)/sum(fxpairraw{i,3:4}); % get rewards per total moves

            elseif str2double(daysess(k).runnum)<7
                for rt=1:2
                    % my last well ~= my current well and my last well ~=
                    % his current well
                    candidates=daysess(k).ratsamples{rt}.thiswell~=daysess(k).ratsamples{rt}.lastwell & ...
                        daysess(k).ratsamples{rt}.hiswell~=daysess(k).ratsamples{rt}.lastwell;
                    % my last well ~= my current well and my current well
                    % == his current well
                    wins=daysess(k).ratsamples{rt}.thiswell~=daysess(k).ratsamples{rt}.lastwell & ...
                        daysess(k).ratsamples{rt}.hiswell==daysess(k).ratsamples{rt}.thiswell;
                    % % of candidates are wins
                    combosraw.(varnames{rt})(i,cumct)=sum(wins)/sum(candidates); % get performance rate
                    combosraw.(varnames{rt+2})(i,cumct)=sum(candidates); % get total moves
                end
                % rews / tot trans
                combosraw.totMatches(i,cumct)=sum(daysess(k).ratsamples{rt}.match)/(combosraw{i,3}(cumct)+combosraw{i,4}(cumct)); % get matches per total moves
                combosraw.totRewards(i,cumct)=sum(daysess(k).ratsamples{rt}.Reward)/(combosraw{i,3}(cumct)+combosraw{i,4}(cumct)); % get rewards per total moves

                cumct=cumct+1;
            end
        end
    end

    
    myvars=combosraw.Properties.VariableNames;
    combomeans=table;
    for i=1:length(myvars)
        combomeans.(myvars{i})=mean(combosraw.(myvars{i}),2);
    end

    
    
    
    subplot(2,3,iterator);
    % arm transitions go to friends arm
    if ~runtots
        % lets fit these data to a binomial probability distribution
        % binofit (wins tots, p0
        [ctrlsraw.binomean, ctrlsraw.binoBounds] = binofit(ctrlsraw.perfA.*ctrlsraw.totA+ctrlsraw.perfB.*ctrlsraw.totB,...
            ctrlsraw.totA+ctrlsraw.totB,0.5);
        [fxpairraw.binomean, fxpairraw.binoBounds] = binofit(fxpairraw.perfA.*fxpairraw.totA+fxpairraw.perfB.*fxpairraw.totB,...
            fxpairraw.totA+fxpairraw.totB,0.5);
        % combos is more complicated
        [combosraw.binomean, combosraw.binoBounds] = binofit(round(combomeans.perfA.*combomeans.totA+combomeans.perfB.*combomeans.totB),...
            round(combomeans.totA+combomeans.totB),0.5);
        
        
        plot(ctrlsraw.binomean); hold on; plot(fxpairraw.binomean); plot(mean(combosraw.binomean,3));
        hold on;
        mycolors=lines(3);
        
        patch([1:size(ctrlsraw,1) size(ctrlsraw,1):-1:1]',[ctrlsraw.binoBounds(:,1)' flipud(ctrlsraw.binoBounds(:,2))'],mycolors(1,:),...
            'FaceAlpha',.5,'EdgeColor','none');
        patch([1:size(fxpairraw,1) size(fxpairraw,1):-1:1]',[fxpairraw.binoBounds(:,1)' flipud(fxpairraw.binoBounds(:,2))'],mycolors(2,:),...
            'FaceAlpha',.5,'EdgeColor','none');
        patch([1:size(combosraw,1) size(combosraw,1):-1:1]',[combosraw.binoBounds(:,1)' flipud(combosraw.binoBounds(:,2))'],mycolors(3,:),...
            'FaceAlpha',.5,'EdgeColor','none');
        hold on;
        plot([7 7],[.2 .9],'k');
        xlabel('Session number');
        ylabel(sprintf('Likelihood of transitioning \n to peer-occupied arm'));
        legend('Ctrl-Ctrl pair','FX-FX pair','FX-Ctrl pair');
        box off
        title(sprintf('cc-fxfx p=%.3e, cc-fxc p=%.3e, \n fxfx-fxc p=%.2e',...
            ranksum(mean(ctrlsraw{:,1:2},2),mean(fxpairraw{:,1:2},2)),...
            ranksum(mean(ctrlsraw{:,1:2},2),mean(combomeans{:,1:2},2)),...
            ranksum(mean(combomeans{:,1:2},2),mean(fxpairraw{:,1:2},2))));
        
        
        barcolors=lines(3);
        % a paired bar graph
        minheight=min([height(ctrlsraw) height(combomeans) height(fxpairraw)]);
        meanvals=[mean(ctrlsraw{1:minheight,1:2},2) mean(combomeans{1:minheight,1:2},2) mean(fxpairraw{1:minheight,1:2},2)]-.5;
        
        
        
    else
        
        
        % and this is total rewards per total arm transitions
        
        % cohort 1 had a wonky day on the 11th day
        if cohort==1 && height(ctrlsraw)==23
            ctrlsraw([11 21 22 23],:)=[]; fxpairraw([11 21 22 23],:)=[]; combomeans([11 21 22 23],:)=[];
        end
        plot((ctrlsraw.totMatches)); hold on;
        plot((fxpairraw.totMatches));
        plot((combomeans.totMatches));
        % lets fit these data to a binomial probability distribution
        % binofit (wins tots, p0
        [ctrlsraw.binomean, ctrlsraw.binoBounds] = binofit(round(ctrlsraw.totMatches.*(ctrlsraw.totA+ctrlsraw.totB)),...
            ctrlsraw.totA+ctrlsraw.totB,0.5);
        [fxpairraw.binomean, fxpairraw.binoBounds] = binofit(round(fxpairraw.totMatches.*(fxpairraw.totA+fxpairraw.totB)),...
            fxpairraw.totA+fxpairraw.totB,0.5);
        % combos is more complicated
        [combomeans.binomean, combomeans.binoBounds] = binofit(round(combomeans.totMatches.*(combomeans.totA+combomeans.totB)),...
            round(combomeans.totA+combomeans.totB),0.5);
        
        
        plot(ctrlsraw.binomean); hold on; plot(fxpairraw.binomean); plot(mean(combomeans.binomean,3));
        hold on;
        mycolors=lines(3);
        
        patch([1:size(ctrlsraw,1) size(ctrlsraw,1):-1:1]',[ctrlsraw.binoBounds(:,1)' flipud(ctrlsraw.binoBounds(:,2))'],mycolors(1,:),...
            'FaceAlpha',.5,'EdgeColor','none');
        patch([1:size(fxpairraw,1) size(fxpairraw,1):-1:1]',[fxpairraw.binoBounds(:,1)' flipud(fxpairraw.binoBounds(:,2))'],mycolors(2,:),...
            'FaceAlpha',.5,'EdgeColor','none');
        patch([1:size(combomeans,1) size(combomeans,1):-1:1]',[combomeans.binoBounds(:,1)' flipud(combomeans.binoBounds(:,2))'],mycolors(3,:),...
            'FaceAlpha',.5,'EdgeColor','none');
        hold on;
        plot([7 7],[.2 .9],'k');
        xlabel('Session number');
        ylabel(sprintf('Total Matched Pokes \n over \n Total Arm Transitions'));
        legend('Ctrl-Ctrl pair','FX-FX pair','FX-Ctrl pair');
        box off
        title(sprintf('cc-fxfx p=%.3e, cc-fxc p=%.3e, \n fxfx-fxc p=%.2e',...
            ranksum(mean(ctrlsraw{:,1:2},2),mean(fxpairraw{:,1:2},2)),...
            ranksum(mean(ctrlsraw{:,1:2},2),mean(combomeans{:,1:2},2)),...
            ranksum(mean(combomeans{:,1:2},2),mean(fxpairraw{:,1:2},2))));
        
        minheight=min([height(ctrlsraw) height(combomeans) height(fxpairraw)]);
        meanvals=[ctrlsraw.totMatches(1:minheight) combomeans.totMatches(1:minheight) fxpairraw.totMatches(1:minheight) ]-.5;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% bar graph now %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    linkaxes(get(gcf,'Children'),'y');
    

    
    
    
    
    
    % for  bar plot on that cohort instead of aggregated
    %{    
    subplot(2,2,iterator+2);
    bp=bar([1 2 3],mean(meanvals),'FaceColor','flat','EdgeColor','none','FaceAlpha',.8); bp.CData=barcolors([1 3 2],:); 
    hold on;
    plot(repmat([1 2 3],size(ctrlsraw,1),1)',meanvals','k');
    errorbar([1 2 3],mean(meanvals),SEM(meanvals,1),'k');
    set(gca,'XAxisLocation','origin','YTick',[-.2:.1:.3],...
        'YTickLabel',cellfun(@(a) num2str(a), mat2cell(.3:.1:.8,1,ones(6,1)),'UniformOutput',false));
    box off;
    set(gca,'XTick',[]);
    if cohort==1, set(gca,'Ylim',[-.25 .35]); else, set(gca,'Ylim',[-.1 .2]); end
    ylabel(sprintf('Likelihood of transitioning \n to peer occupied arm'));
    %}
    
    allcohorts=[allcohorts; meanvals];
    iterator=iterator+1;
end

subplot(2,4,[6 7]);
%{
grps=repmat(1:3,size(allcohorts,1),1);
boxScatterplot(allcohorts(:)+.5,grps,'position',[1])
hold on;
plot(grps',allcohorts'+.5,'color',[.7 .7 .7]);
plot([0.5 3.5],[.5 .5],'r');
set(gca,'XTickLabel',{'Ctrl-Ctrl','FX-Ctrl','FX-FX'});
ylabel(sprintf('Likelihood of transitioning \n to peer-occupied arm'));

%}
% scrub early sessions...
critcohorts=allcohorts; % criterion days
critcohorts([1:9 21:32],:)=[];

bp=bar([1 2 3],mean(critcohorts),'FaceColor','flat','EdgeColor','none','FaceAlpha',.8); bp.CData=mycolors([1 3 2],:); 
hold on;
errorbar([1 2 3],mean(critcohorts),std(critcohorts,1),'k.');
tickvec=[-.2:.05:.3];
set(gca,'XAxisLocation','origin','YTick',tickvec,...
    'YTickLabel',cellfun(@(a) num2str(a), mat2cell(tickvec+.5,1,ones(length(tickvec),1)),...
    'UniformOutput',false));
box off;
set(gca,'XTick',[]);
ylabel(sprintf('Likelihood of transitioning \n to peer occupied arm'));
xlabel('WT-WT     WT-FX     FX-FX');

[a]=ranksum(critcohorts(:,1),critcohorts(:,2));
a(2)=ranksum(critcohorts(:,1),critcohorts(:,3));
a(3)=ranksum(critcohorts(:,2),critcohorts(:,3));
title(sprintf('FX-FX vs FX-Ctrl p=%.2e, FX-FX vs Ctrl-Ctrl p=%.2e, \n FX-Ctrl vs Ctrl-Ctrl p=%.2e',...
    a(3), a(2), a(1)));
%% lets add error bars to those rates with a binomial distribution
 

% this is actually done above
%{
figure;
subplot(1,2,1);
myconf=cell2mat(ctrlsraw(:,3));
errorbar(1:length(myconf), myconf(:,1), myconf(:,2)-myconf(:,1), myconf(:,3)-myconf(:,1))
hold on;
myconf=cell2mat(fxpairraw(:,3));
errorbar(1:length(myconf), myconf(:,1), myconf(:,2)-myconf(:,1), myconf(:,3)-myconf(:,1))

errorbar(1:length(combostat), combostat(:,1), combostat(:,2)-combostat(:,1), combostat(:,3)-combostat(:,1))

xlabel('Session number');
ylabel(sprintf('Proportion of arm transitions \n to peer-occupied arm'));
legend('Ctrl-Ctrl pair','FX-FX pair','FX-Ctrl pair');

subplot(1,2,2);
% this is the bar graph with error bars
[a,b]=binofit(sum(sum(cellfun(@(a) sum(a), ctrlsraw(:,1:2)))),sum(sum(cellfun(@(a) length(a), ctrlsraw(:,1:2)))));
[a(2),b(2,:)]=binofit(sum(sum(cellfun(@(a) sum(a), fxpairraw(:,1:2)))),sum(sum(cellfun(@(a) length(a), fxpairraw(:,1:2)))));
[a(3),b(3,:)]=binofit(sum(linearize(cellfun(@(a) sum(a), combosraw))),sum(linearize(cellfun(@(a) length(a), combosraw))));

br=bar([1 3 2],a);
hold on;
mycolors=lines(3);
br.FaceColor='flat';
br.CData=mycolors([1 3 2],:);
hold on;
errorbar(1,a(1), b(1,1)-a(1), b(1,2)-a(1),'kx');
errorbar(3,a(2), b(2,1)-a(2), b(2,2)-a(2),'kx');
errorbar(2,a(3), b(3,1)-a(3), b(3,2)-a(3),'kx');
set(gca,'XTick',[1 2 3],'XTickLabel',{'Ctrl-Ctrl','FX-Ctrl','FX-FX'});
ylabel(sprintf('Proportion of arm transitions \n to peer-occupied arm'));
xlabel('Animal Pairing')


% now to use a patch instead


for i=1:length(ctrlsraw)
    [ctrlsraw{i,3}(1), ctrlsraw{i,3}(2:3)] = binofit(sum(ctrlsraw{i,1})+sum(ctrlsraw{i,2}),length(ctrlsraw{i,1})+length(ctrlsraw{i,2}));
    [fxpairraw{i,3}(1), fxpairraw{i,3}(2:3)] = binofit(sum(fxpairraw{i,1})+sum(fxpairraw{i,2}),length(fxpairraw{i,1})+length(fxpairraw{i,2}));
    % combos is more complicated
    combotots=sum(linearize(cellfun(@(a) sum(a), combosraw(i,:,:))));
    combotots(2)=sum(linearize(cellfun(@(a) length(a), combosraw(i,:,:))));
    [combostat(i,1), combostat(i,2:3)] = binofit(combotots(1),combotots(2));  
end


figure;
mycolors=lines(3);
subplot(1,2,1);
myconf=cell2mat(ctrlsraw(:,3));
plot(1:length(myconf), myconf(:,1),'Color',mycolors(1,:));
hold on;
hp=patch([1:length(myconf) fliplr(1:length(myconf))]', [myconf(:,2); flipud(myconf(:,3))],...
    mycolors(1,:),'FaceAlpha',.5,'EdgeColor','none');

myconf=cell2mat(fxpairraw(:,3));
plot(1:length(myconf), myconf(:,1),'Color',mycolors(2,:));
hold on;
hp(2)=patch([1:length(myconf) fliplr(1:length(myconf))]', [myconf(:,2); flipud(myconf(:,3))],...
    mycolors(2,:),'FaceAlpha',.5,'EdgeColor','none');


%errorbar(1:length(combostat), combostat(:,1), combostat(:,2)-combostat(:,1), combostat(:,3)-combostat(:,1))
plot(1:length(combostat), combostat(:,1),'Color',mycolors(3,:));
hold on;

hp(3)=patch([1:length(combostat) fliplr(1:length(combostat))]', [combostat(:,2); flipud(combostat(:,3))],...
    mycolors(3,:),'FaceAlpha',.5,'EdgeColor','none');

xlabel('Behavior Day'); box off; xlim([0 21]);
ylabel(sprintf('Proportion of arm transitions \n to peer-occupied arm'));
legend(hp,'Ctrl-Ctrl pair','FX-FX pair','FX-Ctrl pair');
legend('box','off');


subplot(1,2,2);
% this is the bar graph with error bars
[a,b]=binofit(sum(sum(cellfun(@(a) sum(a), ctrlsraw(10:end,1:2)))),sum(sum(cellfun(@(a) length(a), ctrlsraw(10:end,1:2)))));
[a(2),b(2,:)]=binofit(sum(sum(cellfun(@(a) sum(a), fxpairraw(10:end,1:2)))),sum(sum(cellfun(@(a) length(a), fxpairraw(10:end,1:2)))));
[a(3),b(3,:)]=binofit(sum(linearize(cellfun(@(a) sum(a), combosraw(10:end,:,:)))),sum(linearize(cellfun(@(a) length(a), combosraw(10:end,:,:)))));

br=bar([1 3 2],a);
hold on;

br.FaceColor='flat';
br.CData=mycolors([1 3 2],:);
hold on;
errorbar(1,a(1), b(1,1)-a(1), b(1,2)-a(1),'kx');
errorbar(3,a(2), b(2,1)-a(2), b(2,2)-a(2),'kx');
errorbar(2,a(3), b(3,1)-a(3), b(3,2)-a(3),'kx');
set(gca,'XTick',[1 2 3],'XTickLabel',{'Ctrl-Ctrl','FX-Ctrl','FX-FX'});
ylabel(sprintf('Proportion of arm transitions \n to peer-occupied arm'));
xlabel('Animal Pairing')

%}
%%
% i wonder if their armtransitions per minute were higher.
genotypetable=table([1 4; 1 2; 1 2; 2 3],[2 3; 3 4; 3 4; 1 4], 'VariableNames',{'controls','fx'});
genotable={'ctrl','fx','fx','ctrl';...
    'ctrl','fx','ctrl','fx';...
    'ctrl','ctrl','fx','fx';...
    'fx','ctrl','ctrl','fx'};


for i=1:length(ratinfo)
    
    ratinfo(i).genotypes=genotable(ratinfo(i).cohortnum,ratinfo(i).ratnums);
   
end

for k=[1 2 3]
firsthit=[]; fxpair=[];
cumct=[1 1];

    cohortinfo=ratinfo(cell2mat({ratinfo.cohortnum})==k);

for i=1:length(cohortinfo)
    if any(cohortinfo(i).ratnums==1) &&...
            any(cohortinfo(i).ratnums==4) &&...
            str2double(cohortinfo(i).runnum)<7
        % session start, session end
        firsthit(cumct(1),1)=cohortinfo(i).ratsamples{1}.start(1);
        firsthit(cumct(1),2)=cohortinfo(i).ratsamples{1}.end(end);
        % and now the number of transitions that day
        firsthit(cumct(1),3)=sum(diff(cohortinfo(i).ratsamples{1}.thiswell)~=0)+...
            sum(diff(cohortinfo(i).ratsamples{2}.thiswell)~=0);
        cumct(1)=cumct(1)+1;
    elseif any(cohortinfo(i).ratnums==2) &&...
            any(cohortinfo(i).ratnums==3) && ...
            str2double(cohortinfo(i).runnum)<7
        fxpair(cumct(2),1)=cohortinfo(i).ratsamples{1}.start(1);
        fxpair(cumct(2),2)=cohortinfo(i).ratsamples{1}.end(end);
        % and now the number of rewards
        fxpair(cumct(2),3)=sum(diff(cohortinfo(i).ratsamples{1}.thiswell)~=0)+...
            sum(diff(cohortinfo(i).ratsamples{2}.thiswell)~=0);
        cumct(2)=cumct(2)+1;
    end
end

%firsthit([12 17 23],:)=[];
%fxpair([12 17 23],:)=[];
subplot(2,2,k);
plot(firsthit(:,3)./(firsthit(:,2)-firsthit(:,1)).*60); hold on;
plot(fxpair(:,3)./(fxpair(:,2)-fxpair(:,1)).*60);
xlabel('Session number');
ylabel('Number of arm transitions per minute');
legend('WT-WT pair','FX-FX pair');

title(sprintf('ranksum test p= %.2e',...
    ranksum(firsthit(:,3)./max(firsthit(:,1:2),[],2),fxpair(:,3)./max(fxpair(:,1:2),[],2))));
end
%% so one question is what is the arm bias and arm transition bias of these animals?


% do animals visit one arm over the others
% (max(armvisits)-~max(armvisits)/tot visits)

% do animals have a transition bias?

% for i=1:3
%    transbias(i)=diff trans I to a, and trans I to b over sum
% end

for i=1:length(ratinfo)
    % arm bias forst
    for j=1:2
        % armvisits (this well id)
        transitions=accumarray(ratinfo(i).ratsamples{j}.thiswell([false; diff(ratinfo(i).ratsamples{j}.thiswell)~=0]),1);
        ratinfo(i).armbias(j)=max(transitions)/mean(transitions);
        
        % last well id of transitions
        fromtrans=ratinfo(i).ratsamples{j}.lastwell([false; diff(ratinfo(i).ratsamples{j}.thiswell)~=0]);
        totrans=ratinfo(i).ratsamples{j}.thiswell([false; diff(ratinfo(i).ratsamples{j}.thiswell)~=0]);
        for k=1:3
            transitions=[accumarray(totrans(fromtrans==k),1); 0; 0];
            if length(transitions)>3
                ratinfo(i).armtransbias(j,k)=max(transitions(1:3~=k))/mean(transitions(1:3~=k));
            else
               ratinfo(i).armtransbias(j,k)=nan; 
            end
        end
        
        % samples in which I'm going to a new well, and his well is not
        % where I just left
        candidates=ratinfo(i).ratsamples{j}.thiswell~=ratinfo(i).ratsamples{j}.lastwell & ...
            ratinfo(i).ratsamples{j}.hiswell~=ratinfo(i).ratsamples{j}.lastwell;
        % my last well ~= my current well and my current well
        % == his current well
        wins=ratinfo(i).ratsamples{j}.thiswell~=ratinfo(i).ratsamples{j}.lastwell & ...
            ratinfo(i).ratsamples{j}.hiswell==ratinfo(i).ratsamples{j}.thiswell;
        
        ratinfo(i).sessperf(j)=sum(wins)/sum(candidates);
        
    end
end

%% run glm on dwell times 
% with rewarded, whether theres a peer there, whether there isnt...
% or whether you just left a friends well


%% this splits each session into two runs
clear allruns;
allruns.date=ratinfo(1).date;
allruns.datenum=ratinfo(1).datenum;

allruns.ratnum=ratinfo(1).ratnum(1);
allruns.runnum=ratinfo(1).runnum;
allruns.samples=ratinfo(1).ratsamples{1};
allruns(2).date=ratinfo(1).date;
allruns(2).ratnum=ratinfo(1).ratnum(2);
allruns(2).runnum=ratinfo(1).runnum;
allruns(2).samples=ratinfo(1).ratsamples{2};
allruns(2).datenum=ratinfo(1).datenum;

cumct=3;
for i=2:length(ratinfo)
    allruns(cumct).date=ratinfo(i).date;
    allruns(cumct).ratnum=ratinfo(i).ratnum(1);
    allruns(cumct).runnum=ratinfo(i).runnum;
    allruns(cumct).samples=ratinfo(i).ratsamples{1};
    allruns(cumct).datenum=ratinfo(i).datenum;
    cumct=cumct+1;
    allruns(cumct).date=ratinfo(i).date;
    allruns(cumct).ratnum=ratinfo(i).ratnum(2);
    allruns(cumct).runnum=ratinfo(i).runnum;
    allruns(cumct).samples=ratinfo(i).ratsamples{2};
    allruns(cumct).datenum=ratinfo(i).datenum;
    cumct=cumct+1;
end

    
%%  Do any animals show a high bias?
biasnames={}; biasgeno={};
biases=[];
cumct=1;
for i=1:length(ratinfo)
    for j=1:2
        if any(ratinfo(i).cohortnum==[1 3 4])
            biasnames{cumct}=sprintf('%s%d', ratinfo(i).cohortname,ratinfo(i).ratnums(j));
            biasgeno{cumct}=sprintf('%s%d%s', ratinfo(i).cohortname,ratinfo(i).ratnums(j),ratinfo(i).genotypes{j});
            biases(cumct,1)=ratinfo(i).armbias(j);
            biases(cumct,2)=mean(ratinfo(i).armtransbias(j,:));
            cumct=cumct+1;
        end
    end
end

[a,b,c]=anovan(biases(:,1),{biasgeno});     figure;      
multcompare(c); title('arm bias');
 figure;      
[a,b,c]=anovan(biases(:,2),{biasgeno});           
multcompare(c); title('arm transition bias');
    

    
    
    
    
    
    
    