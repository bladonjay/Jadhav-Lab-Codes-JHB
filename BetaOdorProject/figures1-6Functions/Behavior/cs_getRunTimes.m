%Get run times- 
% time between nosepoke off and reward start
% correct trials only
% align nosepoke off times with reward times, take time in between
% nx2 matrix with start and stop

clear
[topDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS44'};

for a= 1:length(animals)
    animal = animals{a};
    animDir = [topDir, animal,'Expt\',animal,'_direct\'];
    
    runeps = cs_getRunEpochs(animDir, animal,'odorplace');
    
    days = unique(runeps(:,1));
    
    odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',days);
     nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',days);
     rewards = loaddatastruct(animDir, animal,'rewards',days);
    
     for day = days'
         daystr = getTwoDigitNumber(day);
         eps = runeps((runeps(:,1) == day),2);
         for ep = eps'
             [cl,cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
             npoff = nosepokeWindow{day}{ep}([cl;cr],2);
             rew = [rewards{day}{ep}.leftWindows(:,1);rewards{day}{ep}.rightWindows(:,1)];
             
             if length(npoff) ~= length(rew)
                 keyboard
             end
             
             run = [npoff,rew];
             runTimes{day}{ep} = run;
             
         end
         
         save([animDir, animal,'runTimes',daystr],'runTimes');
         clear runTimes
     end
end
     