%cs_betaPower_noOdor


close all
clear
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS39','CS41','CS42','CS44'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};
figure
hold on
labels_all = [];
results_all = [];
pvals = [];
bootstrap = 1;
freq = 'resp';
calc = 1;
switch freq
    case 'beta'
        bandpass = [15 30];
    case 'resp'
        bandpass = [7 8];
end

for r = 1:length(regions)
    odorplace = [];
    noOdor = [];

    region = regions{r};
    disp(['Doing ', region])
    if calc ~= 1
        load([topDir,'AnalysesAcrossAnimals\',freq,'Coh_noOdor'])
        odorplace = coh_noOdor.(erase(region,'-'))(coh_noOdor.(erase(region,'-'))(:,2) == 1,1);
        noOdor = coh_noOdor.(erase(region,'-'))(coh_noOdor.(erase(region,'-'))(:,2) == 0,1);
    else
        
    for a = 1:length(animals)
        
        animal = animals{a};
        disp(['Doing animal ',animal]);
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        tetinfo = loaddatastruct(animDir, animal,'tetinfo');
        npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
        odorTriggers = loaddatastruct(animDir, animal,'odorTriggers');
        rewards = loaddatastruct(animDir, animal,'rewards');
        % odor place trials
        daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');

        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = daymatrix(daymatrix(:,1) == day,2);
            
            load([animDir, animal, 'coherence',region,daystr]);
%             tetfilt = ['strcmp($area,''',region,''')'];
%             tets = evaluatefilter(tetinfo,tetfilt);

            %load([animDir, animal, 'specfloor',region,daystr]);
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                time = coherence{day}{ep}.time;
                data = coherence{day}{ep}.Coh;
                
                switch freq
                        case 'beta'
                            
                            data = coherence{day}{ep}.rawCoh;
                            sd = coherence{day}{ep}.rawSD;
                            mn = coherence{day}{ep}.rawMean;
                            data = (data-mn)./sd;
                        case 'resp'
                            
                            data = coherence{day}{ep}.rawCoh;
                            time = coherence{day}{ep}.time;
                            %get baseline for zscore
                            rtimes = [rewards{day}{ep}.leftWindows;rewards{day}{ep}.rightWindows];
                            rewarddata = data(:,isExcluded(time,rtimes));
                            mn = mean(rewarddata,2);
                            sd = std(rewarddata,0,2);
                
                            %zscore to reward coh
                            data = (data-mn)./sd;
                    end
                %freq = coherence{day}{ep}.freq;

            goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
            data = data(find(goodrows),:);
            data = mean(data,1);
         
                
                wins = npWins{day}{ep};
                
                %correct trials only?? 
                [cl, cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                wins = wins([cl;cr],:);
                %get beta power for each trial
                for t = 1:size(wins)
                    d = mean(data(isExcluded(time,wins(t,:))));
                    odorplace = [odorplace;d];
                end
                
            end
            
        end
        
        % no Odor trials
        daymatrix = cs_getRunEpochs(animDir, animal, 'noodor');

        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = daymatrix(daymatrix(:,1) == day,2);
            load([animDir, animal, 'coherence',region,daystr]);

%             tetfilt = ['strcmp($area,''',region,''')'];
%             tets = evaluatefilter(tetinfo,tetfilt);

            %load([animDir, animal, 'specfloor',region,daystr]);
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                
                time = coherence{day}{ep}.time;
                data = coherence{day}{ep}.Coh;
                goodrows = coherence{day}{ep}.freq >= bandpass(1) & coherence{day}{ep}.freq <= bandpass(2);
            data = data(find(goodrows),:);
            data = mean(data,1);
         
                wins = npWins{day}{ep};
                [cl, cr] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                wins = wins([cl;cr],:);
                %get beta power for each trial
                for t = 1:size(wins)
                    d = mean(data(isExcluded(time,wins(t,:))));
                    noOdor = [noOdor;d];
                end
                
            end
            
        end
        %%
        
    end
    
    end
    
    odorplace = odorplace(~isnan(odorplace));
    Cerr = stderr(odorplace);
    
    noOdor = noOdor(~isnan(noOdor));
    Ierr = stderr(noOdor);
    
    if bootstrap == 1
        iterations = 1000;
        Cdist = [];
        for i = 1:iterations
            samp = datasample(1:size(odorplace,1),size(noOdor,1));
            trials = odorplace(samp);
            Cdist = [Cdist;mean(trials)];
        end
        p_upper  = sum(Cdist >= mean(noOdor)) /length(Cdist);
        p_lower = sum(Cdist <= mean(noOdor))/length(Cdist);
        p = min([p_upper, p_lower]);
                pvals = [pvals;p];

        %disp([region, ' pvalue = ', num2str(p)]);
        %text(r-0.1,0.2,['p = ', num2str(p)]);
    CI = getCI(Cdist,95);
        
        plot([r r], CI','k-')
        hold on
        plot(r,mean(Cdist),'ko')
        plot(r + 0.25,mean(noOdor),'ro');

        %xlim([0 3])
        %ylim([-0.1 0.3])
        %title(region);
        %ylabel('beta coherence')
        text(r, mean(noOdor),['p = ',num2str(round(p,2,'significant'))]);
        
        results = [odorplace;noOdor];
        labels = [ones(length(odorplace),1);zeros(length(noOdor),1)];
        coh_noOdor.(erase(region,'-')) = [results,labels];
        box off
                disp([region, ' pvalue = ', num2str(p)]);
                
%         errorbar(r, mean(odorplace), stderr(odorplace),'b.')
%         hold on
%         errorbar(r+0.25, mean(Cdist), 2*std(Cdist),'r.')
        
        %plot(r, mean(C),'k.');
        %axis([0 r+1.25 0.1 0.75])
%         data{r}.C = mean(C);
%         data{r}.I = Idist;
        
%         results = [odorplace;Cdist];
%         labels = [repmat([region(1),'c'],length(odorplace),1); repmat([region(1),'i'],length(Cdist),1)];
%         results_all = [results_all;results];
%         labels_all = [labels_all;labels];
        %plotvars = [plotvars; results,labels];
    
    else 
         %calc significance
        %[~,p] = ttest2(C,I);
        p = ranksum(odorplace,noOdor);
        pvals = [pvals;p];
        disp([region, ' pvalue = ', num2str(p)]);
        
        errorbar(r, mean(odorplace), stderr(odorplace),'b.')
        hold on
        errorbar(r+0.25, mean(noOdor), stderr(noOdor),'r.')
        %plot
%         errorbar(r,mean(C),Cerr,'k.');
%         hold on
%         errorbar(r+0.25,mean(I),Ierr,'k.');
%         data{r}.C = C;
%         data{r}.I = I;
        results = [odorplace;noOdor];
        labels = [repmat([region(1),'c'],length(odorplace),1); repmat([region(1),'i'],length(noOdor),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        
    end
    
    
    %rrPowerData.(region) = C;
    %save([topDir,'AnalysesAcrossAnimals\rrPowerData'],'rrPowerData');
end
save([topDir,'AnalysesAcrossAnimals\',freq,'Coh_noOdor'],'coh_noOdor');
xticks([1.125 2.125 3.125]);
xticklabels(regions);

axis([0.75 3.75 -0.3 0.4]);
% text(1.5, -.2, ['p = ', num2str(round(pvals(1),2))]);
% text(2.5, -.2, ['p = ', num2str(round(pvals(2),2))]);
% text(3.5, -.2, ['p = ', num2str(round(pvals(3),2))]);
ylabel([freq,' Coherence']);
figfile = [figDir, 'EEG\',freq,'Coherence_noOdor'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile); 

% 