
function [rtimes, nptimes, dur] = cs_alignNPReward(rtimes,nptimes,dur)
%all 3 inputs should be 1d vectors

% this function is used to make sure all np times are aligned with the
% appropriate reward. if a np or reward is missing, it is removed from
% analysis.

if rtimes(1) < nptimes(1)
    rtimes(1) = [];
end

maxlength = min(length(nptimes), length(rtimes));
for i = 2:maxlength
    
    
    if nptimes(i) > rtimes(i)
        rtimes(i) = [];
        if length(nptimes) == length(rtimes)
            break
        end
    end
    
    if nptimes(i+1) < rtimes(i) 
        nptimes(i) = [];
        dur(i) = [];
        if length(nptimes) == length(rtimes)
            break
        end
    end
    
    if i == size(nptimes,1) && i < size(rtimes,1)
        rtimes(i+1,:) = [];
    end
    if i == size(rtimes,1) && i < size(nptimes,1)
        nptimes(i+1,:) = [];
        dur(i+1) = [];
    end
    
end