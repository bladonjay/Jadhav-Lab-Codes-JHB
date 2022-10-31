function [eegspecvar] = iseegspecvar(varname)
% eegvar = iseegvar(varname)
% returns 1 if varname is the name of an eeg variable and 0 otherwise
% 
% eeg variables are
% 	'eeg'
%	'theta'
%	'gamma'
%	'ripple'
switch (varname)
case {'eegrefspeclow','eeggndspeclow','eegrefspecmid','eeggndspecmid','eeggndspecfloor','eegrefspechigh'}
    eegspecvar = 1;
otherwise
    eegspecvar = 0;
end

