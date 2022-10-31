%cs_nppos
%plot position, and a dot for each location when rat was in NP
clear
topDir = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS44'};
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt',filesep,animal,'_direct',filesep];
    runepochs = cs_getRunEpochs(animDir, animal, 'odorplace');
        
    days = unique(runepochs(:,1));
    for day = days'
        %load([animDir,animal,'rippletimes']);
        pos = loaddatastruct(animDir, animal,'pos',day);
        nosepoke = loaddatastruct(animDir, animal,'nosepokeWindow',day);
        epochs = runepochs(runepochs(:,1) == day,2);
        for ep = epochs'
        %find pos for each np
        np = nosepoke{day}{ep};
        postime = pos{day}{ep}.data(:,1);
        posdata = pos{day}{ep}.data(:,[2,3]);
        nppos = posdata(isExcluded(postime,np),:);
        
        plot(pos{day}{ep}.data(:,2),pos{day}{ep}.data(:,3));
        hold on
        plot(nppos(:,1),nppos(:,2),'r.')
        
        title([animal,' day ',num2str(day),' epoch ',num2str(ep)])
        pause
        close
        end
        
    end
end