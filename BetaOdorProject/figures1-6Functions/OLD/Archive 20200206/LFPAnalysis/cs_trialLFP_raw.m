%cs_trialLFP
close all
clear
%plots raw LFP, beta filtered LFP, and theta/HRR filtered LFP in one figure
%regions = {'CA1','PFC','OB','TC'};
animals = {'CS41','CS42'};
%trigs = 'nosepokeWindow'
[topDir, figDir] = cs_setPaths();
win = [0.5 1.5];

for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir,animal,'Expt\', animal,'_direct\'];
        
        load([animDir, animal, 'tetinfo.mat']);
        dayfiles = dir([animDir, animal, 'nosepokeWindow*']);
        
        for d = 1:length(dayfiles)
            load(dayfiles(d).name);
            day = length(nosepokeWindow);
            daystr = getTwoDigitNumber(day);
            epochs = find(~cellfun(@(x) isempty(x), nosepokeWindow{day}));
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                epstr = getTwoDigitNumber(epoch);
                windows = nosepokeWindow{day}{epoch};
                exits = windows(:,2);
                windows = [windows(:,1) - win(1), windows(:,1) + win(2)];
                
                %CA1 tetrode
                tetfilter = 'isequal($area,''CA1'')';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                numcells = cellfetch(tetinfo{day}{epoch}(tets),'numcells');
                [~,ind] = max(cell2mat(numcells.values));
                ca1tet = tets(ind(1));
                tetstr = getTwoDigitNumber(ca1tet);
                ca1eeg = loadeegstruct(animDir, animal, 'eeg',day,epoch,ca1tet);
                
                %PFC tetrode
                tetfilter = 'isequal($area,''PFC'')';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                numcells = cellfetch(tetinfo{day}{epoch}(tets),'numcells');
                [~,ind] = max(cell2mat(numcells.values));
                pfctet = tets(ind(1));
                
                pfceeg = loadeegstruct(animDir, animal, 'eeg',day,epoch,pfctet);
                
                
                %OB tetrode
                tetfilter = 'isequal($area,''OB'')';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                obtet = tets(1);
                
                obeeg = loadeegstruct(animDir, animal, 'eeg',day,epoch,obtet);
                obbeta = loadeegstruct(animDir, animal, 'beta',day,epoch,obtet);
                
                %Thermocouple
                tetfilter = 'isequal($area,''TC'')';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                tctet = tets(1);
                
                tceeg = loadeegstruct(animDir, animal, 'eeg',day,epoch,tctet);
                
                times = geteegtimes(obeeg{day}{epoch}{obtet});
                
                
                for tr = 1:size(windows,1)
                    trialinds = isExcluded(times,windows(tr,:));
                    
                    trialtime = times(trialinds);
                    exittime = exits(tr)-trialtime(1)-win(1);
                    trialtime = trialtime - trialtime(1)-win(1);
                    
                    
                    ca1 = ca1eeg{day}{epoch}{ca1tet}.data(trialinds);
                    pfc = pfceeg{day}{epoch}{pfctet}.data(trialinds);
                    ob = obeeg{day}{epoch}{obtet}.data(trialinds);
                    tc = tceeg{day}{epoch}{tctet}.data(trialinds);
                    beta = obbeta{day}{epoch}{obtet}.data(trialinds);
                    
                    figure,
                    set(gcf, 'Position', [1000, 80 ,1000,900]);
                    subplot(5, 1, 1), hold on
                    plot(trialtime, tc, 'k-');
                    plot([0 0], [min(tc(:)) max(tc(:))], 'k--');
                    plot([exittime exittime], [min(tc(:)) max(tc(:))], 'k--');
                    title('Thermocouple')
                    
                    subplot(5, 1, 2), hold on
                    plot(trialtime, ob, 'g-');
                    plot([0 0], [min(ob(:)) max(ob(:))], 'k--');
                    plot([exittime exittime], [min(ob(:)) max(ob(:))], 'k--');
                    title('OB')
                    
                    subplot(5, 1, 3), hold on
                    plot(trialtime, ca1, 'b-');
                    plot([0 0], [min(ca1(:)) max(ca1(:))], 'k--');
                    plot([exittime exittime], [min(ca1(:)) max(ca1(:))], 'k--');
                    title('CA1')
                    
                    subplot(5, 1, 4), hold on
                    plot(trialtime, pfc, 'r-');
                    plot([0 0], [min(pfc(:)) max(pfc(:))], 'k--');
                    plot([exittime exittime], [min(pfc(:)) max(pfc(:))], 'k--');
                    title('PFC')
                    
                    subplot(5, 1, 5), hold on
                    plot(trialtime, beta, 'g-');
                    plot([0 0], [min(beta(:)) max(beta(:))], 'k--');
                    plot([exittime exittime], [min(beta(:)) max(beta(:))], 'k--');
                    title('Beta Filtered')
                    
                    suptitle([animal, ' day ',daystr,' epoch ',epstr,' trial ',num2str(tr)]);
                    pause,
                    close all
                end
                
            end
        end
end