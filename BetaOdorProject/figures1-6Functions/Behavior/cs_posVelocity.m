%cs_posVelocity

%align to lin pos

clear
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS44'};
%win = [1 7];

 %xi = linspace(-win(1),win(2),200); 
 %ti = interp1([-win(1):win(2)],[-win(1):win(2)],xi);
posbins = 0:130;
velocity = [];

%vel_all = [];
for a= 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    %odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',days);
%     nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',days);
    pos = loaddatastruct(animDir, animal, 'pos',days);
    linpos = loaddatastruct(animDir, animal,'linpos',days);
    %velocity = [];
    for r = 1:size(runeps,1)
        ep = runeps(r,:);
        
        %time = linpos{ep(1)}{ep(2)}.statematrix.time;
        traj = linpos{ep(1)}{ep(2)}.statematrix.traj; %only use outbound trials (traj 1 and 3)
        vel = pos{ep(1)}{ep(2)}.data(:,5);
        lindist = linpos{ep(1)}{ep(2)}.statematrix.lindist  ; 
        
        goodinds = (traj == 1 | traj == 3);
        vel = vel(goodinds);
        lindist = lindist(goodinds);
        
            bins = zeros(1,length(posbins-1));
        for b = 1:length(posbins)-1
            binvel = vel(isExcluded(lindist,[posbins(b),posbins(b+1)])); %find velocity of positions that fall within position bin
            
            if isempty(binvel)
                bins(b) = NaN;
            else
            mnbinvel = mean(binvel);
            bins(b) = mnbinvel;
            end
        end
        velocity = [velocity;bins];
        
%         trigs = odorTriggers{ep(1)}{ep(2)}.allTriggers;
%         wins = [trigs-win(1), trigs+win(2)];
        
%         for w = 1:size(wins,1)
%             v = vel(isExcluded(time,wins(w,:)))';
%             
%             %interp to ensure all vectors are the same length
%             xi = linspace(1,length(v),200); 
%             vi = interp1(1:length(v),v,xi);    
%             velocity = [velocity;vi];
%         end
    end
    %vel_all = [vel_all;mean(velocity,1)];
end

mnvel = nanmean(velocity,1);
%err = stderr(velocity);

%  mnvel = mean(vel_all,1);
% err = stderr(vel_all);

figure,
%plot(repmat(ti,5,1)',vel_all');
plot(posbins, mnvel);

axis([0 max(posbins(1:end-1)) 0 40])
xlabel('Linear Distance from NP (cm)');
ylabel('Velocity (cm/second)');

figfile = [figDir,'Behavior\positionVelocity'];
print('-djpeg',figfile);
print('-dpdf',figfile);