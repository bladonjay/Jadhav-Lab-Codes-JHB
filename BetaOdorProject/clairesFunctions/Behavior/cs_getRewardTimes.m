%Creates cell array for each day with reward times.
%Reward Start is when pump is first triggered, 
%Reward End is when reward well is first disengaged after pump turns off. 
clear
animals = {'CS42'};
[topDir,~] = cs_setPaths();

for a = 1: length(animals)
    animal = animals{a};
    dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
    
    cd(dataDir)
    
    taskfiles = dir([animal,'task*']);
    
    epochs = [];
    for tk = 1:length(taskfiles)
        load(taskfiles(tk).name);
        
        epochfilter = '(strcmp($environment, ''odorplace''))';
        eps = evaluatefilter(task,epochfilter);
        epochs = [epochs;eps];
    end
    
    days = unique(epochs(:,1));
    
    for d = 5%1:length(days)
        day = days(d);
        
        if (day<10)
            daystring = ['0',num2str(day)];
        else
            daystring = num2str(day);
        end
        
        runEps = epochs(epochs(:,1) == day, 2);
        
        load([animal,'DIO', daystring, '.mat']);
        load([animal, 'odorTriggers',daystring,'.mat']);
        numEpochs = length(runEps);
        
        for e = 1:length(runEps)
            ep = runEps(e);
            
            if (ep <10)
                epstring = ['0',num2str(ep)];
            else
                epstring = num2str(ep);
            end
            
            %%  Get DIO times and states for reward wells, odor solenoids, nosepoke, and buzzer
            
            [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
            LeftTrigs = odorTriggers{day}{ep}.allTriggers(cl);
            LeftTrigs = [LeftTrigs,ones(length(LeftTrigs),1),nan(length(LeftTrigs),2)];
            
            RightTrigs = odorTriggers{day}{ep}.allTriggers(cr);
            RightTrigs = [RightTrigs,ones(length(RightTrigs),1),nan(length(RightTrigs),2)];
            
            
            
            Din1times = dio{1,day}{1,ep}{1,1}.time;
            Din1states = double(dio{1,day}{1,ep}{1,1}.state);
            leftwellnans = nan([length(Din1states),1]);
            LeftWell = [Din1times, leftwellnans,Din1states, leftwellnans];
            clear('Din1times','Din1hits','leftwellnans');
            
            Din2times = dio{1,day}{1,ep}{1,2}.time;
            Din2states = double(dio{1,day}{1,ep}{1,2}.state);
            rightwellnans = nan([length(Din2states),1]);
            RightWell = [Din2times, rightwellnans,Din2states, rightwellnans];
            clear('Din2times','Din2hits','rightwellnans');
            

            Dout3times = dio{1,day}{1,ep}{1,19}.time;
            Dout3states = double(dio{1,day}{1,ep}{1,19}.state);
            leftpumpnans = nan([length(Dout3states),2]);
            LeftPump = [Dout3times, leftpumpnans, Dout3states];
            clear('Dout3times','Dout3hits','leftpumpnans');
            
            Dout4times = dio{1,day}{1,ep}{1,20}.time;
            Dout4states = double(dio{1,day}{1,ep}{1,20}.state);
            rightpumpnans = nan([length(Dout4states),2]);
            RightPump = [Dout4times, rightpumpnans, Dout4states];
            clear('Dout4times','Dout4hits','rightpumpnans');
            
            
            LeftRewards = sortrows([LeftTrigs; LeftWell; LeftPump]);
            RightRewards = sortrows([RightTrigs; RightWell; RightPump]);
            
            %%

            leftWindows = [];
            rightWindows = [];
            
            leftTrigInds = find(LeftRewards(:,2) == 1);
            for r = 1:length(leftTrigInds)
                
                ind = leftTrigInds(r);
                
                rewardon = find(LeftRewards(ind+1:end,4) == 1,1,'first')+ind;
                
                rewardoff = find(LeftRewards(ind+1:end,4) == 0,1,'first')+ind;
                
                rewardwelloff = find(LeftRewards(rewardoff+1:end,1),1,'first')+rewardoff;
                if isempty(rewardwelloff)
                    rewardwelloff = size(LeftRewards,1);
                end
                
                window = [LeftRewards(rewardon,1) LeftRewards(rewardwelloff, 1)];
                
                leftWindows = [leftWindows; window];
            end
            
            rightTrigInds = find(RightRewards(:,2) == 1);
            for r = 1:length(rightTrigInds)
                
                ind = rightTrigInds(r);
                
                rewardon = find(RightRewards(ind+1:end,4) == 1,1,'first')+ind;
                
                rewardoff = find(RightRewards(ind+1:end,4) == 0,1,'first')+ind;
                
                rewardwelloff = find(RightRewards(rewardoff+1:end,1),1,'first')+rewardoff;
                 if isempty(rewardwelloff)
                    rewardwelloff = size(RightRewards,1);
                 end
                
                window = [RightRewards(rewardon,1) RightRewards(rewardwelloff, 1)];
                
                rightWindows = [rightWindows; window];
            end
            
            windows.leftWindows = leftWindows;
            windows.rightWindows = rightWindows;
            
            clear leftWindows
            clear rightWindows
            
            rewards{1,day}{1,ep} = windows;
            
            
        end
        save([dataDir,animal, 'rewards', daystring],'rewards');
        
        clear rewards
        
    end
    
end
