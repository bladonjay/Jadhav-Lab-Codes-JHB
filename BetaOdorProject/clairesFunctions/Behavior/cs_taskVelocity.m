%cs_taskVelocity

%align to NP DIO time

clear
[topDir figDir] = cs_setPaths();
animals = {'CS44'};
win = [1 7];

 xi = linspace(-win(1),win(2),200); 
 ti = interp1([-win(1):win(2)],[-win(1):win(2)],xi);

velocity = [];
%vel_all = [];
for a= 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',days);
%     nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',days);
    pos = loaddatastruct(animDir, animal,'pos',days);
    %velocity = [];
    v_all = [];
    for r = 1:size(runeps,1)
        ep = runeps(r,:);
        
        time = pos{ep(1)}{ep(2)}.data(:,1);
        vel = pos{ep(1)}{ep(2)}.data(:,5);
        
        trigs = odorTriggers{ep(1)}{ep(2)}.allTriggers;
        wins = [trigs-win(1), trigs+win(2)];
        v_all = [];
        for w = 1:size(wins,1)
            v = vel(isExcluded(time,wins(w,:)))';
            
            %interp to ensure all vectors are the same length
            xi = linspace(1,length(v),200); 
            vi = interp1(1:length(v),v,xi);
            
            v_all= [v_all;vi];
        end
        
    end
    vi_mn = mean(v_all);
        velocity = [velocity;vi_mn];
    %vel_all = [vel_all;mean(velocity,1)];
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
xlabel('Time from odor onset (seconds)');
ylabel('Velocity (cm/second)');

load([topDir, 'AnalysesAcrossAnimals\meanNPduration'])
load([topDir, 'AnalysesAcrossAnimals\rewardLatency'])

plot([meanNP meanNP],[0 max(mnvel)],'k--');
plot([latency latency],[0 max(mnvel)],'k-');

figfile = [figDir,'Behavior\TrialVelocity'];
print('-djpeg',figfile);
print('-dpdf',figfile);