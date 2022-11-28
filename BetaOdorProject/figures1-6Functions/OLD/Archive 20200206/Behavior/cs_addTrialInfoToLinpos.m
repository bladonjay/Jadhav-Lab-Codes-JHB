%Add trialinfo field to linpos- correct/error trials, and r/l odor? 

animals = {'CS31','CS33','CS34','CS35','CS39','CS44'};
topDir = cs_setPaths();

trialInfoKey = "CL = 1, CR = 2, IL = 3, IR =4";

for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    runeps = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(runeps(:,1));
    for d = days'
        odorTriggers = loaddatastruct(animDir,animal,'odorTriggers',d);
        linpos = loaddatastruct(animDir,animal,'linpos',d);
        epochs = runeps(runeps(:,1)==d,2); 
        for ep = epochs'
            [cl, cr, il, ir] = cs_getSpecificTrialTypeInds(odorTriggers{d}{ep});
            time = linpos{d}{ep}.statematrix.time;
            traj = linpos{d}{ep}.statematrix.traj;
            trialinfo = zeros(length(time),1);
            
            %Determine what category each trial falls into, and tag
            %accordingly
            inds = lookup(odorTriggers{d}{ep}.allTriggers,time);
            for i = 1:length(inds)-1
                if ismember(i,cl)
                    trialinfo(inds(i):inds(i+1)-1) = 1;
                elseif ismember(i,cr)
                    trialinfo(inds(i):inds(i+1)-1) = 2;
                elseif ismember(i,il)
                    trialinfo(inds(i):inds(i+1)-1) = 3;
                elseif ismember(i,ir)
                    trialinfo(inds(i):inds(i+1)-1) = 4;
                end
            end
            
            %do last trial separately so that it doesn't error out
            if ismember(inds(end),cl)
                trialinfo(inds(end):end) = 1;
            elseif ismember(inds(end),cr)
                trialinfo(inds(end):end) = 2;
            elseif ismember(inds(end),il)
                trialinfo(inds(end):end) = 3;
            elseif ismember(inds(end),ir)
                trialinfo(inds(end):end) = 4;
            end
            
            
            linpos{d}{ep}.trialinfo = trialinfo;
            linpos{d}{ep}.trialinfokey = trialInfoKey;
        end
        daystr = getTwoDigitNumber(d);
        save([animDir,animal,'linpos',daystr],'linpos');
    end  
end
