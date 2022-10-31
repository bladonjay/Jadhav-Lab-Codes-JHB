%cs_taskVelocity

%align to NP DIO time

%bar graph for pre odor, odor, run, and reward

clear
[topDir figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS44'};
%animals = {'CS34'};
win = [1 7];

 xi = linspace(-win(1),win(2),200); 
 ti = interp1([-win(1):win(2)],[-win(1):win(2)],xi);

v_Odor = [];
v_Pre = [];
v_Run = [];
v_Reward = [];

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
    for r = 2%1:size(runeps,1)
        day = runeps(r,1);
        ep = runeps(r,2);
        
        time = pos{day}{ep}.data(:,1);
        vel = pos{day}{ep}.data(:,5);
        
        [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
        
        %odor 
        wins = nosepokeWindow{day}{ep}([cl;cr],:);
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        v_odor = [v_odor;mean(v)];
        end
        
        %pre-odor 
        diffs = nosepokeWindow{day}{ep}(:,2)-nosepokeWindow{day}{ep}(:,1);
        wins = [nosepokeWindow{day}{ep}([cl;cr],1)-diffs([cl;cr]),nosepokeWindow{day}{ep}([cl;cr],1)];
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        v_pre = [v_pre;mean(v)];
        end
        
        %run
        wins = runTimes{day}{ep};
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        v_run = [v_run;mean(v)];
        end
        
        %reward
        wins = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
        for w = 1:size(wins,1)
        v = vel(isExcluded(time,wins(w,:)))';
        v_reward = [v_reward;mean(v)];
        end
        
  
    end
        v_Pre = [v_Pre;mean(v_pre)];
        v_Odor = [v_Odor;mean(v_odor)];
        v_Run = [v_Run;mean(v_run)];
        v_Reward = [v_Reward;mean(v_reward)];    
end

all = [v_Pre,v_Odor,v_Run,v_Reward];


[p,tbl,stats] = ranova(all);
multcompare(stats)

mnvel = mean(all,1);
err = stderr(all);

figure,
bar(mnvel);
hold on
plot(repmat([1:4],2,1),[mnvel+err;mnvel-err],'k-');

figfile = [figDir,'Behavior\TrialVelocity_bar'];
print('-djpeg',figfile);
print('-dpdf',figfile);