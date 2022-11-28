function out = DFAcs_eventcoherence(index1, index2, excludeperiods, eeg, odorTriggers, varargin)

window = [1.5 1.5];
trigtypes = {'allTriggers'};

movingwin = [1000 20]/1000; 
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 40];
params.tapers = [2 3]; %DO NOT GO LOWER THAN THIS- MESSES UP 

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'trigtypes' 
                trigtypes = varargin{option+1};
            case 'freqband'
                freqband = varargin{option+1}; %i.e. beta, theta
            case 'window'
                window = varargin{option+1}; 
            case 'fpass'
                params.fpass = varargin{option+1};
            
	    otherwise
                error(['Option ',varargin{option},' unknown.']);
        end   
    else
        error('Options must be strings, followed by the variable');
    end
end

odorTriggers = odorTriggers{index1(1)}{index1(2)};

%----- Get the eeg data and time -----%
eeg1 = eeg{index1(1)}{index1(2)}{index1(3)}.data;
time1 = [eeg{index1(1)}{index1(2)}{index1(3)}.starttime:(1/eeg{index1(1)}{index1(2)}{index1(3)}.samprate):eeg{index1(1)}{index1(2)}{index1(3)}.endtime]';
eeg2 = eeg{index2(1)}{index2(2)}{index2(3)}.data;
time2 = [eeg{index2(1)}{index2(2)}{index2(3)}.starttime:(1/eeg{index2(1)}{index2(2)}{index2(3)}.samprate):eeg{index2(1)}{index2(2)}{index2(3)}.endtime]';

        if ~all(time2 == time1)% check if time vectors are matched
              eeg2 = interp1(time2,eeg2,time1,'nearest');
        end

%----- Do the cohereogram calc -----%
disp(['Doing day', num2str(index1(1)), ' epoch', num2str(index1(2)), ' tets ',num2str(index1(3)), '/', num2str(index2(3))]);
[Coh,Phi,~,~,~,t,freq] = cohgramc(eeg1,eeg2,movingwin,params); 
Coh = Coh';
t = t + time1(1);
fs_c = round(1/(t(2)-t(1))); %length of time bins
Nfreq = length(freq); %number of different frequency bands
win = [-window(1):(1/fs_c):window(2)];
winidx = round(win.*fs_c);

%----- Find the mean and std for whole epoch for zscore later -----%
meanCohEpoch = mean(Coh,2);
stdCohEpoch = std(Coh,0,2);
% 
% freqbandindx = (freq >= Hz(1) & freq <= Hz(2))';
% meanCohEpoch_freq = mean(meanCohEpoch(freqbandindx),1); %epoch mean and std within a certain freq range
% stdCohEpoch_freq = stdCohEpoch(freqbandindx);
% 
% 
% mean_meanfreq = mean(meanCohEpoch);
% mean_stdfreq = mean(stdCohEpoch);
% 
% %----- Find cohereogram and coherence in set frequency during trigger
% %window for each trial -----%
% %Does it separately for alltriggers, correcttriggers, and incorrecttiggers 
% 
for g = 1:length(trigtypes)
    trigtype = trigtypes{g};
    triggers = odorTriggers.(trigtype);
    

    allzscore = zeros(length(triggers),length(win));
    allCoh = zeros(Nfreq,length(win));

% 
    for i = 1:length(triggers)
            if triggers(i) > (t(1) + 2) && triggers(i) < (t(end) - 2) % throw away the first  and last 2s
                [junk, trigidx] = min(abs(t - triggers(i)));
                trialindx = [trigidx+winidx];

                %%% COHEROGRAM %%%
                %----- sum coherence within trigger window across trials -----%
                %allCoh = allCoh + Coh(:,trialindx); 
                allCoh(:,:,i) = Coh(:,trialindx);

                %%% COHERENCE LEVEL WITHIN FREQ BAND %%%
                %----- take only the coherence within the frequency band -----%
%                 cohfreqband = Coh(freqbandindx,trialindx);
%                 meanwithinfreqband = mean(cohfreqband,1);

                %----- zscore -----%
%                 meanCohEpoch_rep = repmat(meanCohEpoch_freq,[1, length(cohfreqband)]);
%                 stdCohEpoch_rep = repmat(stdCohEpoch_freq,[1, length(cohfreqband)]);
%                 zscorefreqband = (cohfreqband - meanCohEpoch_rep)./stdCohEpoch_rep;
                
                %meanz = (meanwithinfreqband - mean_meanfreq)./mean_stdfreq;

                %----- mean zscore for the freq band within trigger time window -----% 

                %meanz = mean(zscorefreqband,1);
                %allzscore(i,:) = meanz;

            end
    end

meanAllCoh = mean(allCoh,3);
zscoreCoh = bsxfun(@rdivide,(meanAllCoh - meanCohEpoch),stdCohEpoch(:));


cohgramstr = ['cohgram_',trigtype];
cohfreqstr = [freqband,'coh_',trigtype];

% out.(cohgramstr) = zscoreCoh;
% out.(cohfreqstr) = allzscore;


out.index = [index1,index2(3)];

out.timewindow = win; %saves timewindow. It will save the same thing for every iteration... maybe there is a better way to do this
end