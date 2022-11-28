%sniff rate

%get sniff rate before and after odor

%3/24/20 - changed to do this over sessions rather than trials
clear
[topDir, figDir] = cs_setPaths();
animals = {'CS41','CS42'};

presniffs = [];
postsniffs = [];

presniffs_ses = [];
postsniffs_ses = []; 
for a = 1:length(animals)
    animal = animals{a};
    animDir = [topDir,animal,'Expt\',animal,'_direct\'];
    
    nosepokeWindow = loaddatastruct(animDir,animal,'nosepokeWindow');
    tetinfo = loaddatastruct(animDir,animal,'tetinfo');
    
    filt = 'strcmp($area,''TC'')';
    TCtet = evaluatefilter(tetinfo{1}{1},filt);
    daymat = cs_getRunEpochs(animDir, animal,'odorplace');
    days = unique(daymat(:,1));
    
    for day = days'
        epochs = daymat((daymat(:,1)==day),2);
        
        sniffs_ses = [];
        for ep = epochs'
            windows = nosepokeWindow{day}{ep};
            prewindows = [windows(:,1)-(windows(:,2)-windows(:,1)),windows(:,1)];
            eeg = loadeegstruct(animDir, animal,'eeg',day,ep,TCtet);
            tcdata = eeg{day}{ep}{TCtet}.data(:,1);
            times = geteegtimes(eeg{day}{ep}{TCtet});
           
            %remove times that are too close to end of data
            bad = windows(:,2)>times(end);
           
            windows(bad,:) = [];
            prewindows(bad,:) = [];
            sniffs_ses = [sniffs_ses;zeros(size(windows,1),2)];
            for t = 1:size(windows,1)
                preind = isExcluded(times,prewindows(t,:));
                pre = tcdata(preind);
                p = smoothdata(pre);
                %findpeaks occasionally picks up noise in the data, so downsample
                p = downsample(p,25); 
                s = findpeaks(p);
                sniffs = length(s.loc);
                sniffrate = sniffs/(windows(t,2)-windows(t,1));
                sniffs_ses(t,1) = sniffrate;
                %check to make sure sniff rate isn't impossibly high - if
                %it is, it's probably noise
                if sniffrate > 20
                    figure
                    plot(p)
                    hold on
                    plot(s.loc, p(s.loc),'gs');
                    keyboard
                end
                presniffs = [presniffs;sniffrate];

                postind = isExcluded(times,windows(t,:));
                 post = tcdata(postind);
                p = smoothdata(post);
                p = downsample(p,25);
                s = findpeaks(p);
                sniffs = length(s.loc);
                sniffrate = sniffs/(windows(t,2)-windows(t,1));
                postsniffs = [postsniffs;sniffrate];
                sniffs_ses(t,2) = sniffrate;
                if sniffrate > 20
                    figure
                    plot(p)
                    hold on
                    plot(s.loc, p(s.loc),'gs');
                    keyboard
                end
                
            end
           
        end
         presniffs_ses = [presniffs_ses; mean(sniffs_ses(:,1))];

         postsniffs_ses = [postsniffs_ses; mean(sniffs_ses(:,2))];

    end
end

% p = signrank(presniffs,postsniffs);
% data = [presniffs;postsniffs];
% labels = [zeros(length(presniffs),1); ones(length(postsniffs),1)];
% figure
% boxplot(data,labels);
% text(1,10,['p = ',num2str(p)]);
% xticklabels({'pre','post'});
% ylabel('Sniff Rate');
% ylim([0 12])
% 
p = signrank(presniffs_ses,postsniffs_ses);
data = [presniffs_ses;postsniffs_ses];
labels = [zeros(length(presniffs_ses),1); ones(length(postsniffs_ses),1)];
figure
boxplot(data,labels);
text(1,10,['p = ',num2str(p)]);
xticklabels({'pre','post'});
ylabel('Sniff Rate');
ylim([0 10])
mns = mean([presniffs_ses, postsniffs_ses]);
sems = stderr([presniffs_ses, postsniffs_ses]);

% figure
% errorbar([1 2],[mean(presniffs) mean(postsniffs)],[stderr(presniffs) stderr(postsniffs)])
% ylim([0 8]);
sniffRate.pre = presniffs;
sniffRate.post = postsniffs;

figfile = [figDir,'Behavior\sniffRate'];
%saveas(gcf,newfigfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile);
save([topDir,'AnalysesAcrossAnimals\sniffRate'],'sniffRate')

