%% run_behavior_stats
%
% So i think one way ot analyze these data is to put it into a big GLM and
% run it that way, then we can analyze the impact on the rats behavior both
% at each level and across

% so i dont think this needs to be too complicated, one option is to run a
% glm using loglink and input all the durations and categoricals
% so the table will be...
% 1. study duration
% 2. delay duration (end study to start test sample)
% 3. test1 sample dur
% 4. 
for i=1:length(SuperRat)
    glmTable=array2table(SuperRat(i).Behavior.trialMat)
end

    


%% the other option is an anova with continuous variables

% the crucial questions are as follows:
% 1. what are the important durations and are they consistent?
%   study time, delay time, sample time, choice likelihood (first choice)
% 2. do they all do better than chance?
    % binomial test
% 3. does anyone do worse or better than peers?
    % friedman test one way
% 4. do any durations differ across animals
    % could cat all trials and run by animal, but probably session means
    % could help- this will show that one rat had single camera and
    % durations are different
% 5. do they affect correctness?
    % probably where the GLM would be helpful, by doing a leave one out-
    % this is because durations will be correlated
% 6. what biases do the rats have
% 6a. side biases
    % what %% of sessions does the rat prefer one side? basically all
% 6b. first sample bias
    % what %% of first samples are 'digs'
    % what %% of errors are first sample digs?
% 6c. streakyness
    % error % after error? and error %% after two errors
    % could get clever and run by errors in last 5 trials, 1 to 5 and plot
    % aggregate performance (it will drop)


% preliminary conclusions are that study time doesnt differ correct
% incorrect, but test time does, and rats appear to make most errors when
% they erroneously dig in the first pot.  this suggests rats make
% comparison errors and are not forgetting...


%% 2. do they all do better than chance?
% &2. does anybody do worse than their peers? yes

% how did all the rats do?
performancemat=nan(48,100); % first 100 trials
for i=1:length(SuperRat)
    importing=SuperRat(i).Behavior.trialMat.sample1digcorr;
    keep=SuperRat(i).Behavior.trialMat.keepmat;
    performancemat(i,1:sum(keep))=importing(keep==1,2);
    % run a binomial test on each session
    chanceprob(i)=myBinomTest(sum(importing(keep==1,2)),sum(keep==1),.5,'one');
    Ratinfo(i)=SuperRat(i).Sessinfo;
    fprintf('%s was %.2e \n',Ratinfo(i).name,chanceprob(i));
end
eachrat=nanmean(performancemat,2);
eachrat(:,2)=[Ratinfo(:).ratNumber];
ratnames=cellfun(@(a) a(1:find(a==' ')-1), {Ratinfo.name},'UniformOutput',false)';

[a,b,c]=kruskalwallis(eachrat(:,1),ratnames,'off');

% for a pretty figure of performance
figure;
scatter(eachrat(:,2)*2,eachrat(:,1),8,[.7 .7 .7],'filled','XJitter','rand');
hold on;
errorbar(unique(eachrat(:,2))*2,accumarray(eachrat(:,2),eachrat(:,1),[],@mean),...
    accumarray(eachrat(:,2),eachrat(:,1),[],@std)*2.04,'.','CapSize',0,...
    'LineWidth',2,'Color','k','MarkerSize',12,'Marker','_')
set(gca,'YLim',[.5 1]);
set(gca,'XTick',[2:2:12],'XTickLabel',unique(ratnames));
ylabel('Percent Correct');
xlabel('Rat');

title(sprintf('kruskalwallis test, Chi-Sq(%d) =%.2e p= %.2e\n session max p binom test = %.2e',...
    b{2,3},b{2,5},b{2,6}),max(chanceprob))

%% 3a. do any durations differ across animals, or do they affect behavioral
% performance within any animal? or should i do session...

for i=1:length(SuperRat)
    % get session means instead here
    rawMat=[behavior.samplestart(:,1), behavior.boardup, behavior.sample1,...
            behavior.sample1digcorr];
        myNames={'StudyTime','DelayTime','TestTime','Dig?','Correct'};
        designMat=[designMat; rawMat(:,2)-rawMat(:,1), rawMat(:,3)-rawMat(:,2),...
            rawMat(:,4)-rawMat(:,3), rawMat(:,[5 6])];  
        % tab this for all sessions
        % and run anova for rat
end
%% 3b/c

% i think because these are covariates, lets run it for each rat and run it
% like a glm where ss is just remove that variable.  It'll be a poisson
% link and the regressors will be all the durations and whether they dug or
% didnt on first sample

% we'll test three durations: study duration, overall delay, and test
% duration, then we'll also add in categorical dig on first
%% 3b do the durations affect behavior by rat
eachrat=cellfun(@(a) a.ratNumber, {SuperRat.Sessinfo}); % get tab of all rats
for i=1:max(eachrat)
    % build a supermatrix of a given rat
    ratSess=find(eachrat==i);
    designMat=[];
    for j=1:length(ratSess)
        behavior=SuperRat(ratSess(j)).Behavior.trialMat;
        rawMat=[behavior.samplestart(:,1), behavior.boardup, behavior.sample1,...
            behavior.sample1digcorr];
        myNames={'StudyTime','DelayTime','TestTime','Dig?','Correct'};
        designMat=[designMat; rawMat(:,2)-rawMat(:,1), rawMat(:,3)-rawMat(:,2),...
            rawMat(:,4)-rawMat(:,3), rawMat(:,[5 6])];  
    end
    designTable=array2table(designMat,'VariableNames',myNames);
    [mdl]=fitglm(designTable,'interactions','CategoricalVars',[4],'Distribution','binomial');
    % drop this into a table so we can reconstruct the likelihoods later

end
%% 3c do any durations affect behavior by session
for i=1:length(SuperRat)
    % build a supermatrix of a given rat
    behavior=SuperRat(ratSess(i)).Behavior.trialMat;
    rawMat=[behavior.samplestart(:,1), behavior.boardup, behavior.sample1,...
        behavior.sample1digcorr];
    myNames={'StudyTime','DelayTime','TestTime','Dig?','Correct'};
    designMat=[rawMat(:,2)-rawMat(:,1), rawMat(:,3)-rawMat(:,2),...
        rawMat(:,4)-rawMat(:,3), rawMat(:,[5 6])];
    designTable=array2table(designMat,'VariableNames',myNames);
    [mdl]=fitglm(designTable,'interactions','CategoricalVars',4,'Distribution','binomial');
    % drop this into a table so we can reconstruct the likelihoods later
end

