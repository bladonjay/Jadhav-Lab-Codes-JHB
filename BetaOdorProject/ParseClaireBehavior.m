function [trialdata,usedfiles,oldbeh] = ParseClaireBehavior(taskfile,daynum)
% function [trialdata,usedfiles,oldbeh] = ParseClaireBehavior(taskfiles)
% this function parses behavioral data from files named in the taskfiles.

% i havent quite gotten the file organixzation right, so that is hardcoded,
% but that will come later.  For now this will suffice
% potential input files:
% 1. runTrajBonuds
% 2. runTrialBounds
% 3. nosepokeWindow
% 4. rewardTimes
% 5. rewards
% 6. odorTriggers


% priority list:
%   Get the odor samples period:
%   1. nosepokewindow
%   2. runTrajBounds
%   3. runTrialBounds

%   Get the right;/left and corr/incorr
%   1. odorTriggers

%   get the reward times
%   2. rewardTimes
%   3. rewards


usedfiles={};
oldbeh={};
i=daynum;

% runtrajbouts has start, end lr10, CI10, sniffstart end, rewardstart end,
% stemin and stem out- so this will allow you to get incorrect reward poke
% times
if exist(taskfile{1},'file')==2 % 'runTrajBounds'
    usedfiles=taskfile{1};
    load(taskfile{1});
    runTrajBounds{i}(cellfun(@(a) isempty(a), runTrajBounds{i}))=[];
    oldbeh=runTrajBounds;

    %  now create a substruct
    temptrial=[];
    for k=1:length(runTrajBounds{i})
        temptrial=[temptrial; runTrajBounds{i}{k}.data];
    end
    % now split by column and label (just take the first day
    myfieldnames=runTrajBounds{i}{1}.fields;
    myfieldnames(myfieldnames=='/')=[];
    myfieldnames(myfieldnames=='(' |myfieldnames==')')=[];
    % parse column names
    allspaces=[1 find(myfieldnames==' ') length(myfieldnames)];
    trialdata.(myfieldnames(allspaces(1):allspaces(2)-1))=temptrial(:,1);
    for k=2:length(allspaces)-1
        trialdata.(myfieldnames(allspaces(k)+1:allspaces(k+1)-1))=temptrial(:,k);
    end
    fprintf('used RunTrajBounds File \n');

    % runTrialBounds has start, end, lr10, ci10, sniffstart, sniffend,
    % rewardstart, rewardend, so you get incorrect reward poke times
elseif exist(taskfile{2},'file')==2  % 'runTrialBounds'
    load(taskfile{2});
    usedfiles=taskfile{2};
    runTrialBounds{i}(cellfun(@(a) isempty(a), runTrialBounds{i}))=[];
    oldbeh=runTrialBounds;

    %  now create a substruct
    temptrial=[];
    for k=1:length(runTrialBounds{i})
        temptrial=[temptrial; runTrialBounds{i}{k}.data];
    end
    % now split by column and label
    myfieldnames=runTrialBounds{i}{1}.fields;
    myfieldnames(myfieldnames=='/')=[];
    myfieldnames(myfieldnames=='(' |myfieldnames==')')=[];
    % parse each column name
    allspaces=[1 find(myfieldnames==' ') length(myfieldnames)];
    trialdata.(myfieldnames(allspaces(1):allspaces(2)-1))=temptrial(:,1);
    for k=2:length(allspaces)-1
        trialdata.(myfieldnames(allspaces(k)+1:allspaces(k+1)-1))=temptrial(:,k);
    end
    fprintf('used RunTrialBounds File \n');

    % if you dont have these full datasets, construct it on your own
else
    % this will get your odor on off times and THATS IT
    if exist(taskfile{3},'file')==2 % nosepokeWindow
        load(taskfile{3});
        usedfiles=taskfile{3};
        nosepokeWindow{i}(cellfun(@(a) isempty(a), nosepokeWindow{i}))=[];
        oldbeh=nosepokeWindow;
        % cat all the sessions
        temptrial=[];
        for k=1:length(nosepokeWindow{i})
            temptrial=[temptrial; nosepokeWindow{i}{k}];
        end
        trialdata.sniffstart=temptrial(:,1);
        trialdata.sniffend=temptrial(:,2);
        fprintf('used nosepokeWindow File \n');
    else
        fprintf('cant find nosepokeWIndow, so need to find nose poke off times \n')
    end
    % odorTriggers is better, it has all, triggers, left triggers, right
    % triggers, correct triggers and incorrect triggers
    if exist(taskfile{6},'file')==2 % odorTriggers
        load(taskfile{6});
        usedfiles=[usedfiles;taskfile(6)];
        odorTriggers{i}(cellfun(@(a) isempty(a), odorTriggers{i}))=[];
        oldbeh=[oldbeh; odorTriggers];
        % cat all the sessions
        temptrial=[]; temptrial2=[];
        for k=1:length(odorTriggers{i}) % have to work these up in each block
            ltrigger=odorTriggers{i}{k}.leftTriggers; % start times
            rtrigger=odorTriggers{i}{k}.rightTriggers;
            ltrigger(:,2)=1; rtrigger(:,2)=0;
            triggerdata=sortrows([rtrigger; ltrigger]);

            ctrigger=odorTriggers{i}{k}.correctTriggers; % start times
            ictrigger=odorTriggers{i}{k}.incorrectTriggers;
            ctrigger(:,2)=1; ictrigger(:,2)=0;

            triggerdata2=sortrows([ctrigger; ictrigger]);

            % cat all sessions
            temptrial=[temptrial; triggerdata];
            temptrial2=[temptrial2; triggerdata2];

        end
        % i hope to god this matches up
        trialdata.leftright10=temptrial(:,2);
        trialdata.CorrIncorr10=temptrial2(:,2);

        fprintf('used odorTriggers File \n');
    end
    % reward times will have most of your data
    % this has allTriggers, left Triggers, rightTriggers, correctTriggers
    % and incorrectTriggers
    if exist(taskfile{4},'file')==2 % rewardTimes
        load(taskfile{4});
        usedfiles=[usedfiles;taskfile(4)];
        rewardTimes{i}(cellfun(@(a) isempty(a), rewardTimes{i}))=[];
        oldbeh=[oldbeh; rewardTimes];
        % cat all the sessions
        temptrial=[];
        for k=1:length(rewardTimes{i}) % have to work these up in each block
            triggerdata=rewardTimes{i}{k}.allTriggers; % start times
            [~,b]=intersect(rewardTimes{i}{k}.allTriggers,rewardTimes{i}{k}.leftTriggers);
            triggerdata(b,2)=1; % second col is leftright10
            [~,b]=intersect(rewardTimes{i}{k}.allTriggers,rewardTimes{i}{k}.correctTriggers);
            triggerdata(b,3)=1; % third is CorrIncorr10
            % cat all sessions
            temptrial=[temptrial; triggerdata];
        end
        trialdata.rewardstart=temptrial(:,1);
        trialdata.leftright10=temptrial(:,2);
        trialdata.CorrIncorr10=temptrial(:,3);
        fprintf('used RewardTimes File \n');

        % THIS DOESNT HAVE REWARDEND TIMES!
    else
        fprintf('No reward times files \n');
    end
    % reards are just the on and off times of the rewarded trials, so
    % you'll miss both the unrewarded trials and the correct vs incorrect
    if exist(taskfile{5},'file')==2 % rewards
        load(taskfile{5});
        usedfiles=[usedfiles;taskfile(5)];
        rewards{i}(cellfun(@(a) isempty(a), rewards{i}))=[];
        oldbeh=[oldbeh; rewards];
        % cat all the sessions
        temptrial=[];
        for k=1:length(rewards{i}) % have to work these up in each block
            ltrigger=rewards{i}{k}.leftWindows; % start times
            rtrigger=rewards{i}{k}.rightWindows;
            ltrigger(:,3)=1; rtrigger(:,3)=0;
            triggerdata=sortrows([ltrigger; rtrigger]);
            % cat all sessions
            temptrial=[temptrial; triggerdata];
        end
        trialdata.rewardstart=temptrial(:,1);
        trialdata.rewardend=temptrial(:,2);
        % now need to go and find the most recent trial before each
        % reward to flag as correct
        trialdata.CorrIncorr10=zeros(length(trialdata.sniffstart),1);
        for k=1:length(temptrial)
            thisOdorPer=find(trialdata.sniffend<temptrial(k,1),1,'last');
            trialdata.CorrIncorr10(thisOdorPer)=1;
        end
        fprintf('used rewards File \n');

    end
    if ~isfield(trialdata,'CorrIncorr10')
        fprintf('couldnt find trialdata, this is a PROBLEM!!!!!!!!!\n')
    end
    if isfield(trialdata,'rewarden')
        trialdata.rewardend=trialdata.rewarden;
        trialdata=rmfield(trialdata,'rewarden');
    end
end % now gather any nosepokes or rewards that you can
% fx
end

