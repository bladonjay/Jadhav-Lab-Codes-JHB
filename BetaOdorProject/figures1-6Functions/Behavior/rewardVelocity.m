%cs_checkRewardSpeed
%check to ensure that animal speed increases just after reward
%disnegagement- similar to NP disengagement. 
%Otherwise, change how reward periods are defined

clear
[topDir figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
win = [5 1];

xi = linspace(-win(1),win(2),200); 
 ti = interp1([-win(1):win(2)],[-win(1):win(2)],xi);
 
 velocity = [];
 
 for a= 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    rewards = loaddatastruct(animDir, animal, 'rewards',days);
%     nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',days);
    pos = loaddatastruct(animDir, animal,'pos',days);
    
    v_all = [];
    for r = 1:size(runeps,1)
        ep = runeps(r,:);
        
        time = pos{ep(1)}{ep(2)}.data(:,1);
        vel = pos{ep(1)}{ep(2)}.data(:,5);
        
        trigs = [rewards{ep(1)}{ep(2)}.leftWindows(:,2); rewards{ep(1)}{ep(2)}.rightWindows(:,2)];
        wins = [trigs-win(1), trigs+win(2)];
        v_all = [];
        for w = 1:size(wins,1)
            if wins(w,2) > time(end)
                continue
            end
            v = vel(isExcluded(time,wins(w,:)))';
            
            %interp to ensure all vectors are the same length
            xi = linspace(1,length(v),200); 
            vi = interp1(1:length(v),v,xi);
            
            v_all= [v_all;vi];
        end
        
    end
     velocity = [velocity;v_all];
 end
 mnvel = mean(velocity,1);
CI = getCI(velocity,95);

%  mnvel = mean(vel_all,1);
 err = stderr(velocity);
mnvel = smoothdata(mnvel,'gaussian',10);
err = smoothdata(err,'gaussian',10);
 
figure,
%plot(repmat(ti,5,1)',vel_all');
plot(ti, mnvel);
%patch([ti fliplr(ti)],[CI(1,:), fliplr(CI(2,:))],'k','FaceAlpha',0.3)
patch([ti fliplr(ti)],[mnvel+err, fliplr(mnvel-err)],'b','FaceAlpha',0.3)
hold on
plot([0 0],[0 max(mnvel)],'k--');
xlabel('Time from reward offset (seconds)');
ylabel('Velocity (cm/second)');