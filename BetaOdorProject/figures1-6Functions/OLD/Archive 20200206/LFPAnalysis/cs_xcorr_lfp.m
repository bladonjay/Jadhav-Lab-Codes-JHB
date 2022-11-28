%cs_crosscorr
%Take an average of the entire trace- don't find peaks
clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[topDir, figDir] = cs_setPaths();

%bin = 0.01; % 10 ms
tmax = 0.05; % max lag
maxlags = tmax*1500;
%sw1 = bin*2; % for smoothing corrln.

rhythm = 'beta';
regions = {'CA1','PFC','OB'};
regionpairs = combnk(regions,2);

for r = 1:size(regionpairs,1)
    allpairs = [];
    
    for a = 1:length(animals)
        animal = animals{a};
        
        %gather data
        animDir = [topDir,animal,'Expt\', animal,'_direct\'];
        tetinfo = loaddatastruct(animDir, animal,'tetinfo');
        highbeta = loaddatastruct(animDir, animal,'highBeta');
        
        taskfiles = loaddatastruct(animDir, animal, 'task');
        
        days = 1:length(taskfiles);
        
        for d = days
            day = days(d);
            %daystring = getTwoDigitNumber(day);
            epochs = evaluatefilter(taskfiles{day},'strcmp($environment,''odorplace'')');
            
            for ep = epochs'
                disp(['Doing ',animal, ' day ', num2str(day),' epoch ',num2str(ep)])
                betatet1 = evaluatefilter(tetinfo{day}{ep},['strcmp($area,''',regionpairs{r,1},''') & strcmp($descrip2,''betatet'')']);
                betatet2 = evaluatefilter(tetinfo{day}{ep},['strcmp($area,''',regionpairs{r,2},''') & strcmp($descrip2,''betatet'')']);
                
                if isempty(betatet1) || isempty(betatet2)
                    continue
                end
                beta1 = loadeegstruct(animDir, animal, 'beta',day,ep,betatet1);
                beta2 = loadeegstruct(animDir, animal, 'beta',day,ep,betatet2);
                
                eegtime = geteegtimes(beta1{day}{ep}{betatet1});
                
                
                unique(eegtime(2:end)-eegtime(1:end-1));
                beta1 = beta1{day}{ep}{betatet1}.data(:,1);
                beta2 = beta2{day}{ep}{betatet2}.data(:,1);
                
                trialtimes = highbeta{day}{ep}.OB;
                
                for tr = 1:size(trialtimes,1)
                    inds = isExcluded(eegtime, trialtimes(tr,:));
                    b1 = beta1(inds);
                    b2 = beta2(inds);
                    
                    [xc,lags] = xcorr(b1, b2, maxlags);
                    allpairs = [allpairs;xc'];
                end
            end
        end 
    end
    mn = mean(allpairs,1);
    sem = stderr(allpairs);
    lags = lags/1500;
    figure,
    plot(lags,mn)
    hold on
    %patch([lags, fliplr(lags)], [mn-sem,fliplr(mn)+fliplr(sem)], 'k','FaceAlpha',0.2,'EdgeColor','none')
    
    plot([0 0], [min(mn) max(mn)], 'k--')
    figtitle = [regionpairs{r,1},'-',regionpairs{r,2},' lfp cross-correlation'];
    figfile = [figDir,'EEG\',figtitle];
    %saveas(gcf,figfile,'fig');
    print('-dpdf', figfile);
    print('-djpeg', figfile);
end