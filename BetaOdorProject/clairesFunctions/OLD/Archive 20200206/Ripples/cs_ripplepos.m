%cs_ripplepos
%plot position, and a dot for each location where there was a ripple
%detected
close all
clear
[topDir,figDir] = cs_setPaths();
g = gaussian2(3,(3*3)); %for smoothing
animals = {'CS33'};
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt',filesep,animal,'_direct',filesep];
    runepochs = cs_getRunEpochs(animDir, animal, 'odorplace');
        
    days = unique(runepochs(:,1));
    for day = 3%days'
        load([animDir,animal,'rippletimes']);
        pos = loaddatastruct(animDir, animal,'pos',day);
        
        
        
        epochs = runepochs(runepochs(:,1) == day,2);
        for ep = epochs'
            
            daystr = getTwoDigitNumber(day);
        load([animDir, animal,'nosepokeWindow',daystr]);
        wins = nosepokeWindow{day}{ep};
            
            
        %find pos for each ripple
        rips = ripple{day}{ep}.starttime;
        riptimeinds = lookup(rips,pos{day}{ep}.data(:,1));
        posdata = pos{day}{ep}.data(:,[2 3]);
        postime = pos{day}{ep}.data(:,1);
        timestep = postime(2) - postime(1);
        
        %outbound traj only
        linstate = DFTFsj_getlinstate(animDir,animal,runepochs, 6);        
        excludeState = getExcludePeriods(postime,(linstate{day}{ep}.state == 1 | linstate{day}{ep}.state == 3));
        riptimeinds = riptimeinds(~isExcluded(postime(riptimeinds), excludeState));


        plot(posdata(:,1),posdata(:,2));
        hold on
        plot(posdata(riptimeinds,1),posdata(riptimeinds,2),'r.')
        set(gca, 'YDir','reverse')
        
        
        
        %%
        title([animal,' day ',num2str(day),' epoch ',num2str(ep)])
        %pause
        %close
        figfile = [figDir,'Ripples\ripsPos_raw'];
    
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        
        minx = floor(min(posdata(:,1))) - 1;
        maxx = ceil(max(posdata(:,1))) + 1;
        binx = (minx:maxx);
        miny = floor(min(posdata(:,2))) - 1;
        maxy = ceil(max(posdata(:,2))) + 1;
        biny = (miny:maxy);
      
%         binx = min(posdata(:,1)):max(posdata(:,1));
%         biny = min(posdata(:,2)):max(posdata(:,2));
        
        [occupancy] = hist2(posdata(:,1), posdata(:,2), binx, biny);
        
        [rippos, BX, BY] = hist2(posdata(riptimeinds,1),posdata(riptimeinds,2), binx, biny);
        
        %occnormrippos = rippos./(timestep*occupancy);
        
        smoothedocc = filter2(g,occupancy);
        
        
        smootheddata = filter2(g, rippos);
        
        lowocc = smoothedocc <= 2;
        
         smootheddata(lowocc) = -1;
            figure
            imagesc(smootheddata);
            colormap(bone);
            colorbar
            
            figfile = [figDir,'Ripples\ripsPos_smoothed'];
    
        print('-djpeg', figfile);
        print('-dpdf', figfile);
        %occnormrippos = rippos./occupancy;
                    %output.spikerate(nonzero) = output.spikes(nonzero) ./(timestep* output.occupancy(nonzero) );

        end
        
    end
end