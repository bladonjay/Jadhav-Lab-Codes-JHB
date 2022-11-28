%fix full epoch spec files - change time to reflect actual eeg time. 
clear
topDir = cs_setPaths();
animals = {'CS31'};

for a = 1:length(animals)
        animal = animals{a};
        %disp(['Doing ',animal]);
        animDir = [topDir, animal,'Expt\',animal,'_direct\'];
        files = dir([animDir, animal, 'speclow*']);
        
        for f = 1:length(files)
            load([animDir,files(f).name]);
            day = length(spec);
            daystr = getTwoDigitNumber(day);
            load([animDir, animal, 'coherenceCA1-PFC',daystr]);
            eps = find(~cellfun(@isempty, spec{day}));
            
            for ep = eps
                %epstr = getTwoDigitNumber(ep);
                %eegfiles = dir([animDir,'EEG\',animal,'eeg',daystr,'-',epstr,'-*']);
                %load([animDir,'EEG\',eegfiles(1).name]);
                spectime = spec{day}{ep}.time;
                %eegtime = eeg{day}{ep}{length(eeg{day}{ep})}.starttime;
                %spectime= spectime + eegtime;
                cohtime = coherence{day}{ep}.time;
                
                %check
                if length(spectime) ~= length(cohtime)
                    error('Time vectors are not equal lengths')
                end
                
                spec{day}{ep}.time = cohtime;
            end
            
            save([animDir,files(f).name],'spec')

        end
        
end

