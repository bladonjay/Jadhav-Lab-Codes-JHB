%plot animal's position, +/- 1 second from np exit, and overlay with dots for nosepoke exit time
%indicate that np exit corresponds to the initiation of movement
close all
clear
%% Params
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
timew = 1; %time before and after trigger

for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [topDir, animal,'Expt\',animal, '_direct\'];
    pos = loaddatastruct(animDir, animal,'pos');
    odorOff = loaddatastruct(animDir, animal, 'odorOffset');
    
    %get days
    days = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(days(:,1));
    
    for day = days'
        epochs = cs_getRunEpochs(animDir, animal, 'odorplace',day);
        epochs = epochs(:,2);
        
        %get epochs
        for ep = epochs'
            figure, hold on
            p = pos{day}{ep};
            postime = p.data(:,1);
            xpos = p.data(:,2);
            ypos = p.data(:,3);
            
            offsets = odorOff{day}{ep};
            wins = [offsets, offsets+timew];
            
            for w = 1:size(wins,1)
                win = wins(w,:);
                tind = postime>= win(1) & postime<= win(2);
                x = xpos(tind);
                y = ypos(tind);
                
                trig = offsets(w);
                trigind = lookup(trig,postime);
                trigx = xpos(trigind);
                trigy = ypos(trigind);
                
                plot(x,y,'k-')
                plot(trigx, trigy, 'ro','MarkerFaceColor','Red');
                %pause
            end
            pause
        end
    end
end