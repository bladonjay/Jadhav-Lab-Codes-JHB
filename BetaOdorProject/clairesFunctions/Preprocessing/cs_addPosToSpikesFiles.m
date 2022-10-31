function cs_addPosToSpikesFiles(animID, dataDir, sessionNum)

%This function adds pos data to existing spikes files. 
%For use with MountainSort, after converting from mountainlab to filterframework format

sessionString = getTwoDigitNumber(sessionNum);
pos = loaddatastruct(dataDir, animID, 'pos', sessionNum);


spikes = loaddatastruct(dataDir, animID, 'spikes', sessionNum);
if isempty(spikes)
    %if spikes file does not exist for that day, move on to next day
    return
end

epochs = find(~cellfun(@isempty,spikes{sessionNum}));

filt = 'isequal($tag, ''accepted'') | isequal($tag,''mua'')';


for e = 1:length(epochs)
    epoch = epochs(e);
    eppos = pos{sessionNum}{epoch}.data;
    cells = evaluatefilter(spikes{sessionNum}{epoch},filt);
    posstart = eppos(1,1);
    posend = eppos(end,1);
    
    for c = 1:size(cells,1)
        cell = cells(c,:);
        epspikes = spikes{sessionNum}{epoch}{cell(1)}{cell(2)}.data;
        
        %eliminate spikes that fall outside pos time bounds
        epspikes = epspikes(find(epspikes(:,1)> posstart & epspikes(:,1)<posend),:);
        
        %find spike position (from mcz_createNQSpikesFiles)
        posindex = lookup(epspikes(:,1), eppos(:,1));% mcz
        spikepos= eppos(posindex,2:4);            % mcz
        
        epspikes(:,2:4) = spikepos;
        epspikes(:,7) = posindex;
        spikes{sessionNum}{epoch}{cell(1)}{cell(2)}.data = epspikes; 
    end
end

cd(dataDir)
save([animID,'spikes',sessionString], 'spikes');