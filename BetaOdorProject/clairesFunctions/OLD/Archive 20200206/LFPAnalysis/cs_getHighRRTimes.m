%gets RR during nosepoke. Amplitude is fairly consistent throughout entire
%odor sampling period, so just take the whole period. 

%CORRECT TRIALS ONLY

clear all
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS44'};

regions = 'OB';
threshold = 0.5;

trialtype = 'correct';

switch trialtype
    case 'correct'
        trialstr = '';
    case 'incorrect'
        trialstr = '_incorrect';
end

bufferwin = [0.2 0.2];

for a = 1:length(animals)
    animal = animals{a};
    
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    cd(animDir)
     
    dayepochmatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
    days = unique(dayepochmatrix(:,1));
    nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow');
    odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers');
    
    
    for d = 1:length(days)
        day = days(d);
        daystr = getTwoDigitNumber(day);
        epochs = dayepochmatrix(dayepochmatrix(:,1) == day, 2);
        
        load([animDir, animal,'specfloor',region,daystr]);
        %disp(['Doing ', animal,' Day', daystr]);
        
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
            
            spectime = spec{day}{epoch}.time;
            Spec = spec{day}{epoch}.Spec';
            freq = spec{day}{epoch}.freq;
            
            %find RR freq rows
            goodrows = find(freq >= 7 & freq <=8);
            rr = Spec(goodrows,:);
            
            
            %get np windows
            inds = isExcluded(spectime, [npon npoff]);
            npTimes = spectime(inds); %time vector within np
            npRR = rr(inds); %RR during all NPs
            
            %zscore
            baseline = mean(rr);
            stdev = std(rr);
            zscoreRR = (npRR-baseline)/stdev;
            
            %find windows where RR power is above threshold
            goodinds = zscoreRR >= threshold;
            RRtimes = vec2list(goodinds,npTimes);
            
            if ~isempty(goodinds) 
                % make sure windows are at least 300ms long, ~4 cycles
                goodInds = find(RRtimes(:,2) - RRtimes(:,1) >= 0.3);
                
                RRtimes = RRtimes(goodInds,:);
                
                highRRTimes = RRtimes;
                
            else
                highRRTimes = [];
            end
            
        end
        
        %betaPhase{1,day}{1,epoch} = betaPhaseAll;
        highRR{1,day}{1,epoch} = highRRTimes;
        numtrials = length(npon);
        numRRTimes = size(highRRTimes,1);
        disp([animal, ' day ',daystr,' epoch ',epstr,' has ', num2str(numtrials),' trials and ', num2str(numRRTimes), ' RR times']);
        
    end
    
    
    
    
    
    filename = [animal,'highRR',trialstr];
    save([animDir,filename],'highRR')
    clear highRR
    
end
