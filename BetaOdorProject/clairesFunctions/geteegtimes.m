function etimes = geteegtimes(eegstruct)
% function etimes = geteegtimes(eegstruct)
% 
% returns in etimes the times corresponding to each of the data samples in the
% eeg structure.  
% eegstruct must have a starttime and a samprate field.

etimes = eegstruct.starttime + ((0:(length(eegstruct.data(:,1))-1))*(1/eegstruct.samprate));
