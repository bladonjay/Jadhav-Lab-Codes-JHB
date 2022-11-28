%cs_trialLFP
close all
clear
%plots raw LFP, beta filtered LFP, and theta/HRR filtered LFP in one figure
%regions = {'CA1','PFC','OB','TC'};
animals = {'CS41','CS42'};
%trigs = 'nosepokeWindow'
[topDir, figDir] = cs_setPaths();
win = [0.5 1];

for a = 1%1:length(animals)
        animal = animals{a};
        animDir = [topDir,animal,'Expt\', animal,'_direct\'];
        
        load([animDir, animal, 'tetinfo.mat']);
        dayfiles = dir([animDir, animal, 'nosepokeWindow*']);
        
        for d = 1%3:length(dayfiles)
            load([animDir,dayfiles(d).name]);
            day = length(nosepokeWindow);
            daystr = getTwoDigitNumber(day);
            epochs = find(~cellfun(@(x) isempty(x), nosepokeWindow{day}));
            
            for epoch = 2%epochs'
                %epoch = epochs(ep);
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
                
                ca1beta = loadeegstruct(animDir, animal, 'beta',day,epoch,ca1tet);
                
                %PFC tetrode
                tetfilter = 'isequal($area,''PFC'')';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                numcells = cellfetch(tetinfo{day}{epoch}(tets),'numcells');
                [~,ind] = max(cell2mat(numcells.values));
                pfctet = tets(ind(1));
                
                pfceeg = loadeegstruct(animDir, animal, 'eeg',day,epoch,pfctet);
                
                pfcbeta = loadeegstruct(animDir, animal, 'beta',day,epoch,pfctet);
                
                
                %OB tetrode
                tetfilter = 'isequal($area,''OB'')';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                obtet = tets(1);
                
                obeeg = loadeegstruct(animDir, animal, 'eeg',day,epoch,obtet);
                obbeta = loadeegstruct(animDir, animal, 'beta',day,epoch,obtet);
                obresp = loadeegstruct(animDir, animal, 'resp',day,epoch,obtet);
                
                %Thermocouple
                tetfilter = 'isequal($area,''TC'')';
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                tctet = tets(1);
                
                tceeg = loadeegstruct(animDir, animal, 'eeg',day,epoch,tctet);
                
                times = geteegtimes(obeeg{day}{epoch}{obtet});
                
                
                for tr = 31%59:size(windows,1)
                    trialinds = isExcluded(times,windows(tr,:));
                    
                    trialtime = times(trialinds);
                    exittime = exits(tr)-trialtime(1)-win(1);
                    trialtime = trialtime - trialtime(1)-win(1);
                    
                    
                    ca1_raw = ca1eeg{day}{epoch}{ca1tet}.data(trialinds);
                    pfc_raw = pfceeg{day}{epoch}{pfctet}.data(trialinds);
                    ob_raw = obeeg{day}{epoch}{obtet}.data(trialinds);
                    ca1_beta = ca1beta{day}{epoch}{ca1tet}.data(trialinds);
                    pfc_beta = pfcbeta{day}{epoch}{pfctet}.data(trialinds);
                    tc = tceeg{day}{epoch}{tctet}.data(trialinds);
                    ob_beta = obbeta{day}{epoch}{obtet}.data(trialinds);
                    ob_resp = obresp{day}{epoch}{obtet}.data(trialinds);
                    
                    figure,
                    set(gcf, 'Position', [100, 50 ,1000,900]);
                    subplot(8, 1, 1), hold on
                    plot(trialtime, tc, 'k-');
                    xlim([-win(1) win(2)])
                    plot([0 0], [min(tc(:)) max(tc(:))], 'k--');
                    plot([exittime exittime], [min(tc(:)) max(tc(:))], 'k--');
                    title('Thermocouple')
                    
                    
                    subplot(8, 1, 2), hold on
                    plot(trialtime, ca1_raw, 'b-');
                    xlim([-win(1) win(2)])
                    plot([0 0], [min(ca1_raw(:)) max(ca1_raw(:))], 'k--');
                    plot([exittime exittime], [min(ca1_raw(:)) max(ca1_raw(:))], 'k--');
                    title('CA1')
                    
                    subplot(8, 1, 3), hold on
                    xlim([-win(1) win(2)])
                    plot(trialtime, ca1_beta, 'b-');
                    plot([0 0], [min(ca1_beta(:)) max(ca1_beta(:))], 'k--');
                    plot([exittime exittime], [min(ca1_beta(:)) max(ca1_beta(:))], 'k--');
                    title('Beta Filtered')
                    
                    subplot(8, 1, 4), hold on
                    xlim([-win(1) win(2)])
                    plot(trialtime, pfc_raw, 'r-');
                    plot([0 0], [min(pfc_raw(:)) max(pfc_raw(:))], 'k--');
                    plot([exittime exittime], [min(pfc_raw(:)) max(pfc_raw(:))], 'k--');
                    title('PFC')
                    
                    subplot(8, 1, 5), hold on
                    xlim([-win(1) win(2)])
                    plot(trialtime, pfc_beta, 'r-');
                    plot([0 0], [min(pfc_beta(:)) max(pfc_beta(:))], 'k--');
                    plot([exittime exittime], [min(pfc_beta(:)) max(pfc_beta(:))], 'k--');
                    title('Beta Filtered')
                    
                    subplot(8, 1, 6), hold on
                    xlim([-win(1) win(2)])
                    plot(trialtime, ob_raw, 'g-');
                    plot([0 0], [min(ob_raw(:)) max(ob_raw(:))], 'k--');
                    plot([exittime exittime], [min(ob_raw(:)) max(ob_raw(:))], 'k--');
                    title('OB')
                    
                    subplot(8, 1, 7), hold on
                    xlim([-win(1) win(2)])
                    plot(trialtime, ob_beta, 'g-');
                    plot([0 0], [min(ob_beta(:)) max(ob_beta(:))], 'k--');
                    plot([exittime exittime], [min(ob_beta(:)) max(ob_beta(:))], 'k--');
                    title('Beta Filtered')
                    
                    subplot(8, 1, 8), hold on
                    xlim([-win(1) win(2)])
                    plot(trialtime, ob_resp, 'g-');
                    plot([0 0], [min(ob_resp(:)) max(ob_resp(:))], 'k--');
                    plot([exittime exittime], [min(ob_resp(:)) max(ob_resp(:))], 'k--');
                    title('RR Filtered')
                    
                    suptitle([animal, ' day ',daystr,' epoch ',epstr,' trial ',num2str(tr)]);
                    pause,
                    close all
                end
                
            end
        end
end