%cs_betaPower_noOdor


close all
clear
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%animals = {'CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};
figure
hold on
labels_all = [];
results_all = [];
pvals = [];
bootstrap = 1;
for r = 1:length(regions)
    odorplace = [];
    noOdor = [];

    region = regions{r};
    disp(['Doing ', region])
    for a = 1:length(animals)
        
        animal = animals{a};
        disp(['Doing animal ',animal]);
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        tetinfo = loaddatastruct(animDir, animal,'tetinfo');
        npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
        
        % odor place trials
        daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');

        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = daymatrix(daymatrix(:,1) == day,2);
            
%             tetfilt = ['strcmp($area,''',region,''')'];
%             tets = evaluatefilter(tetinfo,tetfilt);

            %load([animDir, animal, 'specfloor',region,daystr]);
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                tetfilt = ['strcmp($area,''',region,''') && strcmp($descrip2,''betatet'')'];
                tet = evaluatefilter(tetinfo{day}{ep},tetfilt);
                
%                 tet = cs_getMostCellsTet(animal, day, ep, region);
%                 
%             if isempty(tet) %if no cells, just pick the first one
%                 tet = evaluatefilter(tetinfo{day},['strcmp($area''',region,''')']);
%                 tet = tet(1,2);
%             end
            
                beta = loadeegstruct(animDir, animal, 'beta',day,ep,tet);
                time = geteegtimes(beta{day}{ep}{tet});
                data = double(beta{day}{ep}{tet}.data(:,3));
                mn = mean(data);
                sd = std(data);
                %zscore
                data = (data-mn)/sd;
                
                wins = npWins{day}{ep};
                
               
                
                %get beta power for each trial
                for t = 1:size(wins)
                    beta = mean(data(isExcluded(time,wins(t,:))));
                    odorplace = [odorplace;beta];
                end
                
            end
            
        end
        
        % no Odor trials
        daymatrix = cs_getRunEpochs(animDir, animal, 'noodor');

        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = daymatrix(daymatrix(:,1) == day,2);
            
%             tetfilt = ['strcmp($area,''',region,''')'];
%             tets = evaluatefilter(tetinfo,tetfilt);

            %load([animDir, animal, 'specfloor',region,daystr]);
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                tetfilt = ['strcmp($area,''',region,''') && strcmp($descrip2,''betatet'')'];
                tet = evaluatefilter(tetinfo{day}{ep},tetfilt);
                %tet = cs_getMostCellsTet(animal, day, ep, region);
            if isempty(tet)
                pause
                break
            end
            
                beta = loadeegstruct(animDir, animal, 'beta',day,ep,tet);
                time = geteegtimes(beta{day}{ep}{tet});
                data = double(beta{day}{ep}{tet}.data(:,3));
                mn = mean(data);
                sd = std(data);
                %zscore
                data = (data-mn)/sd;
                
                wins = npWins{day}{ep};
                
               
                
                %get beta power for each trial
                for t = 1:size(wins)
                    beta = mean(data(isExcluded(time,wins(t,:))));
                    noOdor = [noOdor;beta];
                end
                
            end
            
        end
        %%
        
    end
    
    
    
    odorplace = odorplace(~isnan(odorplace));
    Cerr = stderr(odorplace);
    
    noOdor = noOdor(~isnan(noOdor));
    Ierr = stderr(noOdor);
    
    if bootstrap == 1
        iterations = 1000;
        Idist = [];
        for i = 1:iterations
            samp = datasample(1:size(noOdor,1),size(odorplace,1));
            trials = noOdor(samp);
            Idist = [Idist;mean(trials)];
        end
        p_upper  = sum(Idist >= mean(odorplace)) /length(Idist);
        p_lower = sum(Idist <= mean(odorplace))/length(Idist);
                pvals = [pvals;min([p_upper, p_lower])];

        %disp([region, ' pvalue = ', num2str(p)]);
        %text(r-0.1,0.2,['p = ', num2str(p)]);

        errorbar(r, mean(odorplace), stderr(odorplace),'b.')
        hold on
        errorbar(r+0.25, mean(Idist), 2*std(Idist),'r.')
        
        %plot(r, mean(C),'k.');
        %axis([0 r+1.25 0.1 0.75])
%         data{r}.C = mean(C);
%         data{r}.I = Idist;
        
        results = [odorplace;Idist];
        labels = [repmat([region(1),'c'],length(odorplace),1); repmat([region(1),'i'],length(Idist),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        %plotvars = [plotvars; results,labels];
    
    else 
         %calc significance
        %[~,p] = ttest2(C,I);
        p = ranksum(odorplace,noOdor);
        pvals = [pvals;p];
        disp([region, ' pvalue = ', num2str(p)]);
        
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

xticks([1.125 2.125 3.125]);
xticklabels(regions);

axis([0.75 3.75 -0.2 1.2]);
text(1.5, -.2, ['p = ', num2str(round(pvals(1),2))]);
text(2.5, -.2, ['p = ', num2str(round(pvals(2),2))]);
text(3.5, -.2, ['p = ', num2str(round(pvals(3),2))]);
ylabel('Beta (15-30 Hz) Power');
figfile = [figDir, 'EEG\BetaPower_noOdor'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile); 

% figure
% hold on
% boxplot(results_all,labels_all,'Colors','br','symbol','k.');
% xticks([1.5 3.5 5.5 7.5]);
% xticklabels(regions);
% text(1.5, 4, ['p = ', num2str(round(pvals(1),2))]);
% text(3.5, 4, ['p = ', num2str(round(pvals(2),2))]);
% text(5.5, 4, ['p = ', num2str(round(pvals(3),2))]);
% %text(7.5, 7, ['p = ', num2str(round(pvals(4),2))]);
% %axis([0 r+1.25 -0.1 1.1])
% %text(1,0.16,['p = ', num2str(p)]);
% %plot([0 r+1.25],[0 0],'k--');
% ylabel('Beta (15-30 Hz) Power');
% %legend({'Correct', 'Incorrect'},'Location','East')
% 
% figfile = [figDir, 'EEG\BetaPower_noOdor_boxplot'];
%             
%     print('-djpeg', figfile);
%     print('-dpdf', figfile);
