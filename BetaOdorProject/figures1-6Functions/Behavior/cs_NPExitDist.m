%plot distance from center (NP) as a function of time from NP exit. 
    % should show that time when position diverges from center is aligned
    % with NP exit
    
%plot animal's position, +/- 1 second from np exit, and overlay with dots for nosepoke exit time
%indicate that np exit corresponds to the initiation of movement
close all
clear
%% Params
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
timew = .75; %time before and after trigger
avgdist = [];
figure(1)
hold on
 allx = [];
for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [topDir, animal,'Expt\',animal, '_direct\'];
    pos = loaddatastruct(animDir, animal,'pos');
    linpos = loaddatastruct(animDir, animal,'linpos');
    odorOff = loaddatastruct(animDir, animal, 'odorOffset');
    npWins = loaddatastruct(animDir, animal,'nosepokeWindow');
    
    %get days
    days = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(days(:,1));
    
   
    for day = days'
        epochs = cs_getRunEpochs(animDir, animal, 'odorplace',day);
        epochs = epochs(:,2);
         
        dayx = [];
        %get epochs
        for ep = epochs'
            close all
            %figure(1), hold on
            p = pos{day}{ep};
            postime = p.data(:,1);
            xpos = p.data(:,2);
            ypos = p.data(:,3);
            
            %get pos during np wins, take avg value to use as np pos
            npwins = npWins{day}{ep};
            npinds = isExcluded(postime,npwins);
            npx = xpos(npinds);
            npy = ypos(npinds);
            npposx = mean(npx);
            npposy = mean(npy);
            %nppos = linpos{day}{ep}.wellSegmentInfo.wellCoord(1);
            
            offsets = odorOff{day}{ep};
            wins = [offsets-timew, offsets+timew];
            
            %manual entry - less reliable
%             fid = figure;
%             plot(xpos,ypos);
%             %input where the nosepoke is on the position traj
%             [nppos,~] = ginput(1);
%             close(fid);
           
            for w = 1:size(wins,1)
                win = wins(w,:);
                tind = postime>= win(1) & postime<= win(2);
                t = postime(tind)-offsets(w);
                if isempty(t)
                    continue
                end
                
                x = xpos(tind);
                %x = abs(x - npposx);
                y = ypos(tind);
                
                dist = [];
                for j = 1:length(x)
                    d = sqrt((x(j)-npposx)^2 + (y(j)-npposy)^2);
                    dist = [dist;d];
                end
                x = dist;
                
                
%                 trig = offsets(w);
%                 trigind = lookup(trig,postime);
%                 trigx = xpos(trigind);
%                 trigx = abs(trigx - nppos);
%                 trigy = ypos(trigind);
                
                i = [-timew:1/30:timew];
                x_new = interp1(t,x,i);
                
                dayx = [dayx;x_new];
%                 figure(1)
%                 plot(t,x,'k-')
%                 plot([0 0], [0 10], 'r-');
                
             
                %pause
            end
             
%                 title([animal,' day ',num2str(day),' ep ',num2str(ep)])
                %disp('pause')
                
            
           
            %pause
        end
        allx = [allx; nanmean(dayx)];
        

    end
end

avgdist = nanmean(allx);

 figure(2)
 plot(i(2:end-1),allx(:,2:end-1),'k-','Color',[0.7 0.7 0.7]);
 hold on
 plot(i(2:end-1),avgdist(:,2:end-1),'k-','LineWidth',3);
 xlim([i(2) i(end-1)])
 plot([0 0],[0 15],'r--')
 xlabel('Time from odor offset')
ylabel('Distance from NP (cm)')
figfile = [figDir,'Behavior\NPExitDistance'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
