%cs_rippleRate_NPvsReward

%Compare ripple rate between NP period and reward period
%Use position info, rather than DIO. (Can't use CS41/CS42 in this case
%though)

% 4cm/second speed threshold
% dio times for reward wells, find first time after reward trigger when
% velocity > 4, this is cutoff time for reward. Can use something similar
% for NP
clear
[topDir,figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
speedthresh = 4;
rateNP_session = [];
rateNP_trial= [];
rateReward_session = [];
rateReward_trial = [];

for a = 1:length(animals)
        animal = animals{a};
        
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
       
        runepochs = cs_getRunEpochs(animDir, animal, 'odorplace');
        
        days = unique(runepochs(:,1));
         for day = days'
            odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers', day);
            rewards = loaddatastruct(animDir, animal, 'rewards', day);
            %ripple = loaddatastruct(animDir, animal, 'rippletimes');
            load([animDir,animal,'rippletimes']);
            %cellinfo = loaddatastruct(animDir, animal, 'cellinfo');
            pos = loaddatastruct(animDir, animal,'pos',day);

            epochs = runepochs(runepochs(:,1) == day,2);
            
            for ep = epochs'
                index = [day,ep];
                if strcmp(animal,'CS41') && ismember(index,[4,2],'rows')
                    continue
                end
                %find times around NP and ripples - first time where
                %velocity < 4 --> trig --> first time > 4
                %NP times
                
                trigs = odorTriggers{day}{ep}.correctTriggers;
                windows = [];
                speed = pos{day}{ep}.data(:,5);
                postime = pos{day}{ep}.data(:,1);
                for tr = 1:length(trigs)
                    trig = trigs(tr);
                    posind = lookup(trig,postime);
                    %sometimes, rat's head is still moving when NP starts.
                    %In this case, take first time speed drops below
                    %threshold as the start time. 
                    if speed(posind) > speedthresh
                        posind = find(speed(posind:end) <= speedthresh,1,'first')+(posind-1);
                    end
                    start= postime(find(speed(1:posind)>=speedthresh,1,'last'));
                    stop = postime(find(speed(posind:end)>=speedthresh,1,'first')+(posind-1));
                    
                    %occasionally, rat's speed will not drop below 4 during
                    %np. this causes extraction of very long windows
                    %spanning multiple trials. exclude windows that are
                    %longer than 3 seconds long (generous time for full
                    %stop, np, and turn).
                    if ~isempty(start) && ~isempty(stop) && (stop-start) < 3 
                        windows = [windows;start,stop];
                        riptimes = isExcluded(ripple{day}{ep}.starttime,[start,stop]);
                    count = sum(riptimes);
                    rateNP_trial = [rateNP_trial;count];
                    end
                    
                end
                totaltime = sum(windows(:,2)-windows(:,1));
                %find ripples that occured during windows
                riptimes = isExcluded(ripple{day}{ep}.starttime,windows);
                count = sum(riptimes);
                dur = ripple{day}{ep}.endtime(riptimes)-ripple{day}{ep}.starttime(riptimes);
                
                rateNP_session = [rateNP_session;count];
                
                %Reward times
                trigs = [rewards{day}{ep}.rightWindows(:,1);rewards{day}{ep}.leftWindows(:,1)];
                windows = [];
                for tr = 1:length(trigs)
                    trig = trigs(tr);
                    posind = lookup(trig,pos{day}{ep}.data(:,1));
                    if speed(posind) > speedthresh
                        posind = find(speed(posind:end) <= speedthresh,1,'first')+(posind-1);
                    end
                    %Make sure reward time is longer than 1 second
                    start= postime(find(speed(1:posind)>=speedthresh,1,'last'));
                    stop = postime(find(speed(posind:end)>=speedthresh,1,'first')+(posind-1)); 
                    while stop-start < 1
                        %prevposind = posind;
                        posind = find(speed(posind+1:end) <= speedthresh,1,'first')+(posind);
                        start= postime(posind);%postime(find(speed(prevposind:posind)>=speedthresh,1,'last')+(prevposind-1));
                        stop = postime(find(speed(posind:end)>=speedthresh,1,'first')+(posind-1));
                    end
                        
                    if ~isempty(start) && ~isempty(stop)
                        windows = [windows;start,stop];
                        riptimes = isExcluded(ripple{day}{ep}.starttime,[start,stop]);
                    count = sum(riptimes);
                    rateReward_trial = [rateReward_trial;count];
                    end
                    
                end
                totaltime = sum(windows(:,2)-windows(:,1));
                %find ripples that occured during windows
                riptimes = isExcluded(ripple{day}{ep}.starttime,windows);
                count = sum(riptimes);
                dur = ripple{day}{ep}.endtime(riptimes)-ripple{day}{ep}.starttime(riptimes);
                
                rateReward_session = [rateReward_session;count];
         
            end
         end
end
figure
[p] = ranksum(rateNP_session,rateReward_session);
bar([1,2],[mean(rateNP_session),mean(rateReward_session)]);
disp(['session p value = ',num2str(p)])

figure
[p] = ranksum(rateNP_trial,rateReward_trial);
bar([1,2],[mean(rateNP_trial),mean(rateReward_trial)]);
disp(['trial p value = ',num2str(p)])
hold on
errorbar([mean(rateNP_trial),mean(rateReward_trial)],[stderr(rateNP_trial),stderr(rateReward_trial)],'k.')
text(0.75,1.5,['p = ',num2str(p)]);
xticklabels({'Nosepoke','Reward'})
ylabel('Mean ripples per trial')

figfile = [figDir,'Ripples\NPvsRewardRips'];
print('-djpeg', figfile);
print('-dpdf', figfile);
    