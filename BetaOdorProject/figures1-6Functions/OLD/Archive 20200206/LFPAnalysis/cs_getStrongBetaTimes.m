%for each animal, go to direct folder. load odortriggers file. load tetinfo
%file. get tetrode numbers for each region. for each day and epoch (from
%odortriggers), load all beta files for tetrodes for each region, keep the envelope (column 3). during
%time just before the odor trigger (.5 s?) calc the mean and sd beta power
%across all tetrodes. then, after odor trigger, find the moment when the
%avg beta power is greater than 3sd away from the baseline. save this as
%the start time, then find the moment after that when it is back down to
%within 3 sd of baseline , save as end time. save in cell array for the
%day, with start and end times (one for each odor trigger) for each epoch. 

clear all
%topDir = 'D:\OdorPlaceAssociation\';
topDir = 'F:\Data\OdorPlaceAssociation\';
animals = {'CS31','CS33','CS34','CS35'};
regions = {'CA1','PFC','OB'};
trigtypes = {'allTriggers','correctTriggers','incorrectTriggers'};

for a = 1:length(animals)
    animal = animals{a};
    
    dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
    cd(dataDir)
    
    %get tetinfo file
    tetinfofile = dir([animal,'tetinfo.mat']);
    load(tetinfofile.name);
    
    trigfiles = dir([animal,'odorTriggers*']);
    trigfiles = {trigfiles.name};
    
    for f = 1:length(trigfiles)
        load(trigfiles{f});
        day = length(odorTriggers);
        daystr = getTwoDigitNumber(day);
        epochs = find(~cellfun(@isempty,odorTriggers{1,day}));
        
        disp(['Doing ', animal,' Day', daystr]); 
       
       for tr = 1:length(trigtypes)
           trigtype = trigtypes{tr};
           
       for e = 1:length(epochs)
            epoch = epochs(e);
            triggers = odorTriggers{1,day}{1,epoch}.(trigtype);
        
            for r = 1:length(regions)
               region = regions{r};
               tetfilter = ['((strcmp($area, ''',region,''')) && (strcmp($descrip2, ''betatet'')))'];
               tet = evaluatefilter(tetinfo{1,day}{1,epoch}, tetfilter);
               epstr = getTwoDigitNumber(epoch);
               
               
               
               if ~isempty(tet)
               
                   tetstr = getTwoDigitNumber(tet);
                   load([dataDir,'EEG\',animal,'beta',daystr,'-',epstr,'-',tetstr,'.mat']);
                   beta = beta{day}{epoch}{tet};
                   betatimes = [beta.starttime:1/beta.filtersamprate:beta.endtime]';
                   betaphase = double(beta.data(:,2)); %phase
                   betaamp = double(beta.data(:,3)); %amplitude

               
               
               
                   highbeta_region = [];
                   for tr = 1:length(triggers)
                       trigger = triggers(tr);
                       baselinetimes = [(trigger - 0.5), trigger];
                       baselinebeta = betaamp(betatimes > baselinetimes(1) & betatimes < baselinetimes(2));
                       meanBaseline = mean(baselinebeta);
                       stdBaseline = std(baselinebeta);

                       odortime = [trigger, (trigger + 2)];
                       trigbeta = betaamp(betatimes > odortime(1) & betatimes < odortime(2),:);
                       %trigbetaphase = betaphaseAllTets(betatimes > odortime(1) & betatimes < odortime(2),:);
                       trigbetatimes = betatimes(betatimes > odortime(1) & betatimes < odortime(2));

                       meanbeta = mean(trigbeta,2);
                       meanbeta = filter(ones(1,250)/250, 1, meanbeta); %avg across 250 datapoints

                       highbetastartind = find(meanbeta > (meanBaseline+2*(stdBaseline)),1);
                       highbetastart = trigbetatimes(highbetastartind);
                       remaining = meanbeta(highbetastartind:end);
                       remainingtime = trigbetatimes(highbetastartind:end);
                       highbetaendind = find(remaining <= (meanBaseline+2*(stdBaseline)),1);
                       highbetaend = remainingtime(highbetaendind);
                       %endind = find(trigbetatimes == highbetaend);


                       if highbetaend-highbetastart >= 0.1 %high beta times must be at least 100ms long (about 2 full cycles)
                           if ~isempty(highbetastart) && ~isempty(highbetaend)
                                highbeta_region = [highbeta_region; highbetastart, highbetaend];
                           end
                       end
                  
                       
                   
%                    figure, hold on
%                    plot([1:750],baselinebeta)
%                    plot([751:(750+3000)],meanbeta)
                   
                   
                   end
                   highBetaAll.(region) = highbeta_region;
                   betaPhaseAll.(region) = [betatimes, betaphase];
               else
                   highBetaAll.(region) = [];
                   betaPhaseAll.(region) = [];
               end
               
            end
            
            betaPhase{1,day}{1,epoch} = betaPhaseAll;
            highBeta{1,day}{1,epoch}.(trigtype) = highBetaAll;
            
       end
       
       
       
    end
    end
    
    filename = [animal,'highBeta'];
    save([dataDir,filename],'highBeta')
    clear highBeta
    
    filename = [animal,'betaPhase'];
    save([dataDir,filename],'betaPhase')
    clear betaPhase
end