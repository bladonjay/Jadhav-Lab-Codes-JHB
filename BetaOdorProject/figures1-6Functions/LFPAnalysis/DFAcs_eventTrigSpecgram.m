function [out] = DFAcs_eventTrigSpecgram(index, excludeperiods, eeg, eegspec, Triggers, varargin)

win = [1.5 1.5];

params.Fs = 1500;
params.tapers = [1 1]; % Should I put this in or let it use default tapers
params.err = [2 0.05];
params.fpass = [50 100];
movingwin = [1000 20]/1000; cwin = movingwin;

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
 %change the second number in matrix to change step size (changes
 %resolution)
if params.fpass(2) == 300
    savetag = 'high';
    movingwin = [100 10]/1000; cwin = movingwin;
end
if params.fpass(2) == 100
    savetag = 'mid';
    movingwin = [400 20]/1000; cwin = movingwin;
end
if params.fpass(2) == 40
    savetag = 'low';
    %movingwin = [1000 50]/1000; 
    cwin = movingwin;
end
if params.fpass(2) <= 20
    savetag = 'floor';
    %movingwin = [2000 20]/1000; 
    cwin = movingwin;
end

% if params.fpass(1) == 50 && params.fpass(2) == 100
%     savetag = 'high';
%     %movingwin = [2000 20]/1000; 
%     cwin = movingwin;
% end


% SET DATA
% -------------------------------------------

if exist('eegref')
    eeg = eegref{index(1)}{index(2)}{index(3)};
else
    eeg = eeg{index(1)}{index(2)}{index(3)};
end

eegdata = eeg.data;
eegtimes = (eeg.starttime:(1/eeg.samprate):eeg.endtime)';
starttime = eeg.starttime;
endtime = eeg.endtime;


meandayspec = eegspec{index(1)}{1}{index(3)}.meandayspec; % Stored in first epoch
stddayspec = eegspec{index(1)}{1}{index(3)}.stddayspec;

% if strcmp(savetag, 'high')
%     meandayspec = meandayspec(36:end);
%     stddayspec = stddayspec(36:end);
% end


for tt = 1:length(trigtypes)
    trigtype = trigtypes{tt};
    
    triggers = Triggers{index(1)}{index(2)};
    if isstruct(triggers)
        triggers = Triggers{index(1)}{index(2)}.(trigtype);
    end
    

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
    [S, Stime, Sfreq] = mtspecgramtrigc(eegdata,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);
                    Stime = Stime - win(1); 
    
                    
                    %Get the Zscore
                    S = bsxfun(@minus,S,meandayspec); % Can use meandayspec and stddayspec instead
                    S = bsxfun(@rdivide,S,stddayspec);
                    
                    Smean = mean(S,3);
                    
                    
   % Smoothing along freq axis 
%     for i=1:size(Smean,2)
%     winst = i-smwin/2; winend = i+smwin/2;
%     if winst<1, winst=1; end
%     if winend>size(Smean,2), winend = size(Smean,2); end
%     Smean(:,i) = mean(Smean(:,winst:winend),2);
%     end
% 
% % Smoothing along time axis
%     smwin=smwin*2;
%     for i=1:size(Smean,1)
%     winst = i-smwin/2; winend = i+smwin/2;
%     if winst<1, winst=1; end
%     if winend>size(Smean,1), winend = size(Smean,1); end
%     Smean(i,:) = mean(Smean(winst:winend,:),1);
%     
%     end

    
    out.(trigtype).Smean = Smean;
    out.(trigtype).index = index;

end