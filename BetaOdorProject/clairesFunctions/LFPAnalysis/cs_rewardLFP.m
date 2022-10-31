%cs_rewardLFP
close all
clear
%plots raw LFP, beta filtered LFP, and theta/HRR filtered LFP in one figure
regions = {'CA1'};
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42'};

[topDir, figDir] = cs_setPaths();

for r = 1:length(regions)
    region = regions{r};
    for a = 2:length(animals)
        animal = animals{a};
        animDir = [topDir,animal,'Expt\', animal,'_direct\'];
        
        load([animDir, animal, 'tetinfo.mat']);
        dayfiles = dir([animDir, animal, 'rewards*']);
        
        for d = 1:length(dayfiles)
            load(dayfiles(d).name);
            day = length(rewards);
            daystr = getTwoDigitNumber(day);
            epochs = find(~cellfun(@(x) isempty(x), rewards{day}));
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                epstr = getTwoDigitNumber(epoch);
                windows = [rewards{day}{epoch}.leftWindows;rewards{day}{epoch}.rightWindows];
                tetfilter = ['(isequal($area,''',region,''')) && (strcmp($descrip,''riptet''))'];
                tets = evaluatefilter(tetinfo{day}{epoch},tetfilter);
                
                tetstr = getTwoDigitNumber(tets(1));
                load([animDir,'\EEG\', animal, 'eeg',daystr, '-',epstr,'-',tetstr,'.mat']);
                %timeinds = find(any(timestamps >= windows(:,1) & timestamps < windows(:,2),1));
                
                for t = 1%length(tets):-1:1
                    tet = tets(t);
                    tetstr = getTwoDigitNumber(tet);
                    load([animDir,'\EEG\', animal, 'eeg',daystr, '-',epstr,'-',tetstr,'.mat']);
                    timestamps = eeg{day}{epoch}{tet}.timerange(1):1/eeg{day}{epoch}{tet}.samprate:eeg{day}{epoch}{tet}.timerange(2);
                    eeg_new(t,:) = double(eeg{day}{epoch}{tet}.data(:,1));
                    load([animDir,'\EEG\', animal, 'ripple',daystr, '-',epstr,'-',tetstr,'.mat']);
                    ripple_new(t,:) = double(ripple{day}{epoch}{tet}.data(:,1));
                    load([animDir,'\EEG\', animal, 'theta',daystr, '-',epstr,'-',tetstr,'.mat']);
                    theta_new(t,:) = double(theta{day}{epoch}{tet}.data(:,1));
                end
                
                for tr = 1:size(windows,1)
                    trialinds = find(timestamps >= windows(tr,1)-0.2 & timestamps < windows(tr,1)+8);
                    trialtime = timestamps(trialinds) - windows(tr,1);
                    trialeeg = eeg_new(:,trialinds);
                    trialripple = ripple_new(:,trialinds);
                    trialtheta = theta_new(:,trialinds);
                    
                    figure,
                    set(gcf, 'Position', [800,350,700,500]);
                    subplot(3, 1, 1), hold on
                    plot(trialtime, trialeeg, 'k-');
                    plot([0 0], [min(trialeeg(:)) max(trialeeg(:))], 'k--');
                    
                    subplot(3, 1, 2), hold on
                    plot(trialtime, trialripple, 'b-');
                    plot([0 0], [min(trialripple(:)) max(trialripple(:))], 'k--');
                    
                    subplot(3, 1, 3), hold on
                    plot(trialtime, trialtheta, 'r-');
                    plot([0 0], [min(trialtheta(:)) max(trialtheta(:))], 'k--');
                    
                    pause,
                    close all
                end
               clear eeg_new
               clear ripple_new
               clear theta_new
                
            end
        end
    end
end