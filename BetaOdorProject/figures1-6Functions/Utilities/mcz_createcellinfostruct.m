function mcz_createcellinfostruct(animdirect,fileprefix,append)
% createcellinfostruct(animdirect,fileprefix)
% createcellinfostruct(animdirect,fileprefix,append)
%
% This function created a cellinfo file in the animal's data directory.
% For each cell, the spikewidth, mean rate, and number of spikes is saved.  If a cellinfo file
% exists and new data is being added, set APPEND to 1.
%
% MCZ ADD- If day unclustered, make NaNs

if (animdirect(end) ~= filesep)
    animdirect = [animdirect filesep];
end

cellinfo = [];
if (nargin < 3)
    append = 0;
end
if (append)
    try
        load([animdirect, fileprefix,'cellinfo']);
    end
end
spikefiles = dir([animdirect, fileprefix, 'spikes*']);
for i = 1:length(spikefiles)
    load([animdirect, spikefiles(i).name]);
    timerange = cellfetch(spikes, 'timerange');
    for j = 1:size(timerange.index,1)
        d = timerange.index(j,1);
        e = timerange.index(j,2);
        t = timerange.index(j,3);
        c = timerange.index(j,4);
        cellinfo{d}{e}{t}{c}.spikewidth = timerange.values{j};
        try
            spikewidth = spikes{d}{e}{t}{c}.spikewidth;
        catch
            spikewidth = NaN;
        end
        try
            [csi propbursts] = computecsi(spikes{d}{e}{t}{c}.data(:,1), ...
                spikes{d}{e}{t}{c}.data(:,6), 10);
        catch
            csi = NaN;
            propbursts = NaN;
        end

        try
        epochtime = diff(spikes{d}{e}{t}{c}.timerange);
        numspikes = size(spikes{d}{e}{t}{c}.data,1);
        %d, e, t, c
        tag = spikes{d}{e}{t}{c}.tag;
        cellinfo{d}{e}{t}{c}.meanrate = numspikes/epochtime;
        cellinfo{d}{e}{t}{c}.numspikes = numspikes;
        cellinfo{d}{e}{t}{c}.spikewidth = spikewidth;
        cellinfo{d}{e}{t}{c}.csi = csi;
        cellinfo{d}{e}{t}{c}.propbursts = propbursts;
        cellinfo{d}{e}{t}{c}.tag = tag;
        catch
        tag = spikes{d}{e}{t}{c}.tag;
        cellinfo{d}{e}{t}{c}.meanrate = NaN;
        cellinfo{d}{e}{t}{c}.numspikes = NaN;
        cellinfo{d}{e}{t}{c}.spikewidth = NaN;
        cellinfo{d}{e}{t}{c}.csi = NaN;
        cellinfo{d}{e}{t}{c}.propbursts = NaN;
        cellinfo{d}{e}{t}{c}.tag = tag;
        end
    end
end

save([animdirect,fileprefix,'cellinfo'],'cellinfo');
