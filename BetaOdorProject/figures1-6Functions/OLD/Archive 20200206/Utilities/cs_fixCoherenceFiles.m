function cs_fixCoherenceFiles(topDir, animals, regions)
%CA1-PFC coherence files have wrong times in them- replace wih times from eeg files

for a = 1:length(animals)
    animal = animals{a};
    
    coherencefiles = dir([topDir,animal,'Expt\',animal,'_direct\',animal,'coherence',regions,'*']);
    load([topDir,animal,'Expt\',animal,'_direct\',animal,'tetinfo.mat'])
    
    for d = 1:length(coherencefiles) 
            load([topDir,animal,'Expt\',animal,'_direct\',coherencefiles(d).name])
            coherence_temp = coherence{1,d};
            epochs = find(~cellfun(@isempty, coherence_temp));
            
            for ep = 1:length(epochs)
                
                epoch = epochs(ep);
                tets = find(~cellfun(@isempty, tetinfo{d}{epoch}));
                tet = tets(1); %use first existing tetrode for time
                
                
                load([topDir,animal,'Expt\',animal,'_direct\EEG\', animal, 'eeg',getTwoDigitNumber(d),'-',getTwoDigitNumber(epoch),'-',getTwoDigitNumber(tet),'.mat'])
            
                eegtime = [eeg{d}{epoch}{tet}.starttime:(1/eeg{d}{epoch}{tet}.samprate):eeg{d}{epoch}{tet}.endtime]';
                N = length(eegtime);
                starttime = eegtime(1);
                
                time = starttime + (0.5:0.02:(N-749)/1500); %this calculation is taken from cohgramc.m, exactly what is done when making cohgram files
                
                cohdim = size(coherence_temp{epoch}.coherence,2);
                
                %check, pause if not equal
                if cohdim ~= length(time)
                    pause;
                end
                
                coherence_temp{epoch}.time = time;
                coherence{1,d} = coherence_temp;
            end
             
            save([topDir,animal,'Expt\',animal,'_direct\',animal, 'coherence',regions, getTwoDigitNumber(d),'.mat'],'coherence');
    end
end