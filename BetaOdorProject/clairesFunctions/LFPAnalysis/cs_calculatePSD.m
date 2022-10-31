function [out] = cs_calculatePSD(index, excludeperiods, eeg, eegspec, Triggers, varargin)

% win = [1.5 1.5];
% 
 params.Fs = 1500;
% params.tapers = [1 1]; % Should I put this in or let it use default tapers
% params.err = [2 0.05];
% params.fpass = [50 100];
% movingwin = [1000 20]/1000; cwin = movingwin;
inflection = 0;
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
        case 'inflection'
            inflection = varargin{option+1};
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
        triggeredTimeWindowStart = trigtime - max(win);
        triggeredTimeWindowEnd = trigtime + max(win);
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
goodTrials = triggers;
    triggers = goodTrials;
    triggers = triggers-starttime; endtime = endtime - starttime;
    

    %Remove triggering events that are too close to the beginning or end
    while triggers(1)<win(1)
        triggers(1) = [];
    end
    rem = find(triggers + win(2) > endtime);
    
    triggers(rem) = [];
    
            
    
    meandayspec = eegspec{index(1)}{1}{index(3)}.meandayspec; % Stored in first epoch
    stddayspec = eegspec{index(1)}{1}{index(3)}.stddayspec;


    disp(['Doing PSD calculation. Ntriggers = ',num2str(length(triggers)), ' day ',num2str(index(1)),' epoch ',num2str(index(2)),' tet ',num2str(index(3))]);
    [S, t, ~] = mtspecgramtrigc(eegdata,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);
    %dim 1 = time, dim 2 = freq, dim 3 = trials
    
    %Zscore
    S_norm = bsxfun(@minus,S,meandayspec); % Can use meandayspec and stddayspec instead
    S_norm = bsxfun(@rdivide,S,stddayspec);
                    
                    Smean = mean(S,3);
                    Smean_norm = mean(S_norm,3);
                    
    timewin = sum(win);
    times = [-win(1): timewin/size(S,1): win(2)-(timewin/size(S,1))];
    
    
    
    preinds = find(times<inflection);
    postinds = find(times>=inflection);
    
    
    %save both raw and zscored data
    psd_pre = mean(Smean(preinds,:),1);
    psd_post = mean(Smean(postinds,:),1);
    
    
    psd_all_pre = squeeze(mean(S(preinds,:,:),1));
    psd_all_post = squeeze(mean(S(postinds,:,:),1));

    
    psd_pre_norm = mean(Smean_norm(preinds,:),1);
    psd_post_norm = mean(Smean_norm(postinds,:),1);
    
    psd_all_pre_norm = squeeze(mean(S_norm(preinds,:,:),1));
    psd_all_post_norm = squeeze(mean(S_norm(postinds,:,:),1));
    
    
    %save in output structure
    out.psd_pre = psd_pre;
    out.psd_post = psd_post;
    
    out.psd_all_pre = psd_all_pre;
    out.psd_all_post = psd_all_post;
    
    out.psd_pre_norm = psd_pre_norm;
    out.psd_post_norm = psd_post_norm;
    
    out.psd_all_pre_norm = psd_all_pre_norm;
    out.psd_all_post_norm = psd_all_post_norm;
    
    out.index = index;
