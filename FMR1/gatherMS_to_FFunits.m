%% parsing the spikes data matfiles coming from mountainsort to ff

[spikeFile,spikeDir]=uigetfile('F:\SocialData\Neural\XFB3','Select a spikes.mat file');

load(fullfile(spikeDir,spikeFile));
% you will now have spikes in your workspace

% the data are organized this way
% spikes- 1 by n recording days cell matrix
%    in that its 


sessnum=4;
tetnums=[];
for i=1:length(spikes{sessnum})
    % this is each epoch, so first lets get which cells exist
    tetnums=find([tetnums cellfun(@(a) ~isempty(a), spikes{sessnum}{i})]);
end

% now for each tet, lets gather some info
unitdata=struct;
for i=1:length(tetnums) % for each tetrode
    for j=1:length(spikes{sessnum}) % for each epoch
        tempcell=spikes{sessnum}{j}{tetnums(i)};
        for k=1:length(tempcell)
            if ~isempty(tempcell)
                tempstruct(