function [k_dist, z_dist] = cs_phaselockBootstrap(s, lfptime, phase, wins, nwins_incorr, iterations)

%inputs: 
% s = all spikes for cell
% lfptime = lfptime for the day
% wins = correct windows
% nwins_incorr = number of incorrect trials
% iterations = number of bootsrap iterations

k_dist = [];
z_dist = [];

for i = 1:iterations
    samp = datasample(1:size(wins,1),nwins_incorr);
    trials = wins(samp,:);
    
    sph = phase(lookup(s, lfptime));
    goodspikes = isExcluded(s, trials);
    sph = double(sph(find(goodspikes))) / 10000;  % If no spikes, this will be empty
    if length(sph)> 10 %don't bother with stats unless there are at least 10 spikes to use
        
        % Von Mises Distribution - From Circular Stats toolbox
        [~, k] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
        [~, z] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
        
    else
        k=0;
        z=0;
        
    end
    
    k_dist = [k_dist;k];
    z_dist = [z_dist;z];
    
end
end