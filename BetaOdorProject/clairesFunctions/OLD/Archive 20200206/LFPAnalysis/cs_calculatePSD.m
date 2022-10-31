function [out] = cs_calculatePSD(index, excludeperiods, eeg, Triggers, varargin)

% win = [1.5 1.5];
% 
 params.Fs = 1500;
% params.tapers = [1 1]; % Should I put this in or let it use default tapers
% params.err = [2 0.05];
% params.fpass = [50 100];
% movingwin = [1000 20]/1000; cwin = movingwin;

varargin = varargin{1,1}; %CS 1/8/2017 added so that varargin can be defined once in DFS_eventTrigSpecgram.m and then passed in as a single cell array. Comment out if not using this format.
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'movingwin'
            movingwin = varargin{option+1};
        case 'fpass'
            params.fpass = varargin{option+1};
        case 'trigtypes'
            trigtypes = varargin{option+1};
        case 'gnd'
            do_wrtgnd = varargin{option+1};
        case 'win' 
            win = varargin{option+1};
        case 'tapers'
            params.tapers = varargin{option+1};
    end
end

cwin = movingwin;

% SET DATA
% -------------------------------------------


eeg = eeg{index(1)}{index(2)}{index(3)};
eegdata = eeg.data;
eegtimes = geteegtimes(eeg);
starttime = eeg.starttime;
endtime = eeg.endtime;

    
triggers = Triggers{index(1)}{index(2)}.allTriggers;
    
    

    goodTrials = [];
    for t = 1:length(triggers)  %From cs_removeNoisyTrials
        trigtime = triggers(t,1);
        triggeredTimeWindowStart = trigtime - win(1,1);
        triggeredTimeWindowEnd = trigtime + win(1,2);
        triggeredEEGtimes = eegtimes(eegtimes>triggeredTimeWindowStart & eegtimes<triggeredTimeWindowEnd);
        triggeredEEG = eegdata(eegtimes>triggeredTimeWindowStart & eegtimes<triggeredTimeWindowEnd);

        last = triggeredEEG(2:end);
        first = triggeredEEG(1:end-1);
        diff = last - first;

        %Finds times where there is a large jump between adjacent
        %eeg times - (probably noise)
        bad = find(abs(diff) > 1500); %can change this threshold if necessary

        if isempty(bad) %if there are no large jumps, keep this trial
            goodTrials(end+1,1) = trigtime;
        
        
        end
    end
    triggers = goodTrials;
    triggers = triggers-starttime; endtime = endtime - starttime;
    

    %Remove triggering events that are too close to the beginning or end
    while triggers(1)<win(1)
        triggers(1) = [];
    end
    rem = find(triggers + win(2) > endtime);
    
    triggers(rem) = [];
    
            
    disp(['Doing event-triggered specgram. Ntriggers = ',num2str(length(triggers)), ' day ',num2str(index(1)),' epoch ',num2str(index(2)),' tet ',num2str(index(3))]);
    [S, ~, ~] = mtspecgramtrigc(eegdata,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);

    Smean = mean(S,3);
    psd = mean(Smean,1);
    

    out.psd = psd;
    out.index = index;
