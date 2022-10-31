%for each animal, go to direct folder. load odortriggers file. load tetinfo
%file. get tetrode numbers for each region. for each day and epoch (from
%odortriggers), load all beta files for tetrodes for each region, keep the envelope (column 3). during
%time just before the odor trigger (.5 s?) calc the mean and sd beta power
%across all tetrodes. then, after odor trigger, find the moment when the
%avg beta power is greater than 3sd away from the baseline. save this as
%the start time, then find the moment after that when it is back down to
%within 3 sd of baseline , save as end time. save in cell array for the
%day, with start and end times (one for each odor trigger) for each epoch.


%Define strong beta times differently for each region. Based on average
%spectrograms. E.g.  much stronger beta bursts in OB than in CA1 so
%criteria for high beta in CA1 should be lower?

%CORRECT TRIALS ONLY

clear all
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

% Use OB high beta times only
%regions = {'CA1','PFC','OB'};
region = 'OB';
bufferwin = [0 0];

trialtypes = {'correct','incorrect'};
for tt = 1:length(trialtypes)
    trialtype = trialtypes{tt};
    switch trialtype
        case 'correct'
            trialstr = '';
        case 'incorrect'
            trialstr = '_incorrect';
    end
    
    
    
    for a = 1:length(animals)
        animal = animals{a};
        
        dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
        cd(dataDir)
        
        %get tetinfo file
        tetinfofile = dir([animal,'tetinfo.mat']);
        load(tetinfofile.name);
        
        dayepochmatrix = cs_getRunEpochs(dataDir, animal, 'odorplace');
        days = unique(dayepochmatrix(:,1));
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            epochs = dayepochmatrix(dayepochmatrix(:,1) == day, 2);
            
            load([dataDir, animal, 'nosepokeWindow',daystr,'.mat'])
            load([dataDir, animal, 'odorTriggers',daystr,'.mat'])
            disp(['Doing ', animal,' Day', daystr]);
            
            for e = 1:length(epochs)
                epoch = epochs(e);
                npon = nosepokeWindow{1,day}{1,epoch}(:,1) - bufferwin(1); %add a window around np since beta can start early
                npoff = nosepokeWindow{1,day}{1,epoch}(:,2) + bufferwin(2);
                
                [cl, cr, il, ir] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                switch trialtype
                    case 'correct'
                        trialinds = sort([cl;cr]);
                    case 'incorrect'
                        trialinds = sort([il;ir]);
                end
                npon = npon(trialinds);
                npoff= npoff(trialinds);
                
                epstr = getTwoDigitNumber(epoch);
                
               tet = cs_getMostCellsTet(animal,day,epoch,region);
               tetstr = getTwoDigitNumber(tet);
%                 %for r = 1:length(regions)
%                 %region = regions{r};
%                 tetfilter = ['((strcmp($area, ''',region,''')) && (strcmp($descrip2, ''betatet'')))'];
%                 tet = evaluatefilter(tetinfo{1,day}{1,epoch}, tetfilter);
%                 
                if ~isempty(tet)
                    
                    tetstr = getTwoDigitNumber(tet);
                    load([dataDir,'EEG\',animal,'beta',daystr,'-',epstr,'-',tetstr,'.mat']);
                    beta = beta{day}{epoch}{tet};
                    betatimes = [beta.starttime:1/beta.filtersamprate:beta.endtime]';
                    %betaphase = double(beta.data(:,2)); %phase
                    betaamp = double(beta.data(:,3)); %amplitude
                    
                    baseline = mean(betaamp);
                    stdev = std(betaamp);
                    
                    switch region
                        case 'CA1'
                            threshold = 0.75; %
                        case 'PFC'
                            threshold = 1;
                        case 'OB'
                            threshold = 1;
                    end
                    
                    %find all times during NP
                    betainds = isExcluded(betatimes, [npon, npoff]);
                    npBeta = betaamp(betainds);
                    npTimes = betatimes(betainds);
                    
                    %find beta times that were above threshold
                    zscoreBeta = (npBeta-baseline)/stdev;
                    inds = find(zscoreBeta >= threshold);
                    
                    %separate into distinct periods
                    if ~isempty(inds) %make sure there are some high beta times
                        highBetaTimes = vec2list((zscoreBeta >= threshold), npTimes);
                        % CS- cleaned up 10/31/2019. Was giving some
                        % periods that were too long.
                        
                        %                         endInds = [inds(diff(inds)>1); inds(end)];
                        %                         startInds = [1; endInds(1:end-1)+1];
                        %                         highBetaTimes = [npTimes(startInds), npTimes(endInds)];
                        %                         betaPhase = [npTimes(startInds:endInds), npPhase(startInds:endInds)];
                        
                        % make sure windows are at least 100ms long (~2 beta cycles).
                        goodInds = find(highBetaTimes(:,2) - highBetaTimes(:,1) >= 0.1);
                        
                        highBetaTimes = highBetaTimes(goodInds,:);
                        
                        
                        %end
                        highbeta.(region) = highBetaTimes;
                        %betaPhaseAll.(region) = [highBetaTimes, betaphase];
                    else
                        highbeta.(region) = [];
                    end
                else
                    highbeta.(region) = [];
                    %betaPhaseAll.(region) = [];
                end
                %betaPhase{1,day}{1,epoch} = betaPhaseAll;
                highBeta{1,day}{1,epoch} = highbeta;
                
            end
            
            
            
        end
        
        
        
        filename = [animal,'highBeta',trialstr];
        save([dataDir,filename],'highBeta')
        clear highBeta
        
    end
end

%     filename = [animal,'betaPhase'];
%     save([dataDir,filename],'betaPhase')
%     clear betaPhase