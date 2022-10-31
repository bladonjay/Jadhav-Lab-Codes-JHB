%cs_taskVelocity_avgChange

%align to NP DIO time

%bar graph for pre odor, odor, run, and reward

clear
[topDir figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS44'};
%animals = {'CS34'};
win = [1 7];

 xi = linspace(-win(1),win(2),200); 
 ti = interp1([-win(1):win(2)],[-win(1):win(2)],xi);

% v_Odor = [];
% v_Pre = [];
% v_Run = [];
% v_Reward = [];
% 

preToOdor = [];
odorToRun = [];
runToReward = [];
for a= 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',days);
    nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',days);
    rewards = loaddatastruct(animDir, animal,'rewards',days);
    runTimes = loaddatastruct(animDir, animal,'runTimes',days);
    pos = loaddatastruct(animDir, animal,'pos',days);
    
    v_odor = [];
    v_pre = [];
    v_run = [];
    v_reward = [];
    for r = 1:size(runeps,1)
        day = runeps(r,1);
        ep = runeps(r,2);
        
        time = pos{day}{ep}.data(:,1);
        vel = pos{day}{ep}.data(:,5);
        
        [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
        
        %odor 
        wins = nosepokeWindow{day}{ep}([cl;cr],:);
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        %v_odor = [v_odor;mean(v)];
        v_odor = [v_odor;min(v)];
        end
        
        %pre-odor 
        diffs = nosepokeWindow{day}{ep}(:,2)-nosepokeWindow{day}{ep}(:,1);
        wins = [nosepokeWindow{day}{ep}([cl;cr],1)-diffs([cl;cr]),nosepokeWindow{day}{ep}([cl;cr],1)];
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        %v_pre = [v_pre;v];
        v_pre = [v_pre;max(v)];
        end
        
        %run
        wins = runTimes{day}{ep};
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        %v_run = [v_run;v];
        v_run = [v_run;max(v)];
        end
        
        %reward
        wins = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        %v_reward = [v_reward;v];
        v_reward = [v_reward;min(v)];
        end
        
  
    end
        preToOdor = [preToOdor;mean(v_odor-v_pre)];
        odorToRun = [odorToRun;mean(v_run-v_odor)];
        runToReward = [runToReward;mean(v_reward-v_run)];
        
end
mnchange = mean([preToOdor,odorToRun,runToReward]);
err = stderr([preToOdor,odorToRun,runToReward]);
figure,
bar(mnchange);
hold on
plot(repmat([1:3],2,1),[mnchange+err;mnchange-err],'k-');
xticklabels({'PreToOdor','OdorToRun','RunToReward'})
ylabel('Change in speed');


figfile = [figDir,'Behavior\TrialVelocity_change'];
print('-djpeg',figfile);
print('-dpdf',figfile);