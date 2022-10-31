%cs_betaPowerPerformance

%compare beta power between correct and incorrect trials
%average beta power across nosepoke window, get distribution for each trial
%type

%can also downsample, since there are way fewer incorrect trials
close all
clear
[topDir,figDir] = cs_setPaths();
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};
figure
hold on
labels_all = [];
results_all = [];
pvals = [];
bootstrap = 1;
for r = 1:length(regions)
    C = [];
    I = [];

    region = regions{r};
    disp(['Doing ', region])
    for a = 1:length(animals)
        animal = animals{a};
        disp(['Doing animal ',animal]);
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');
        tetinfo = loaddatastruct(animDir, animal,'tetinfo');
        odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers');
        npWins = loaddatastruct(animDir, animal, 'nosepokeWindow');
        
        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            epochs = daymatrix(daymatrix(:,1) == day,2);
            
%             tetfilt = ['strcmp($area,''',region,''')'];
%             tets = evaluatefilter(tetinfo,tetfilt);

            %load([animDir, animal, 'specfloor',region,daystr]);
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                tet = cs_getMostCellsTet(animal, day, ep, region);
            if isempty(tet)
                break
            end
            
                resp = loadeegstruct(animDir, animal, 'resp',day,ep,tet);
                time = geteegtimes(resp{day}{ep}{tet});
                data = double(resp{day}{ep}{tet}.data(:,3));
                mn = mean(data);
                sd = std(data);
                %zscore
                data = (data-mn)/sd;
                
                [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                correctwins = npWins{day}{ep}([correct_left;correct_right],:);
                incorrectwins = npWins{day}{ep}([incorrect_left;incorrect_right],:);
                
                %get beta power for each trial
                for t = 1:size(correctwins)
                    correctresp = mean(data(isExcluded(time,correctwins(t,:))));
                    C = [C;correctresp];
                end
                for t = 1:size(incorrectwins)
                    incorrectresp = mean(data(isExcluded(time,incorrectwins(t,:))));
                    I = [I;incorrectresp];
                end
            end
            
        end
        
    end
    C = C(~isnan(C));
    Cerr = stderr(C);
    
    I = I(~isnan(I));
    Ierr = stderr(I);
    
    if bootstrap == 1
        iterations = 1000;
        Idist = [];
        for i = 1:iterations
            samp = datasample(1:size(I,1),size(C,1));
            trials = I(samp);
            Idist = [Idist;mean(trials)];
        end
         p_upper  = sum(Idist >= mean(C)) /length(Idist);
        p_lower = sum(Idist <= mean(C))/length(Idist);
                pvals = [pvals;min([p_upper, p_lower])];

        %disp([region, ' pvalue = ', num2str(pvals(end))]);
        %text(r-0.1,0.2,['p = ', num2str(p)]);

        errorbar(r, mean(C), stderr(C),'b.')
        hold on
        errorbar(r+0.25, mean(Idist), 2*std(Idist),'r.')
        
        %plot(r, mean(C),'k.');
        %axis([0 r+1.25 0.1 0.75])
%         data{r}.C = mean(C);
%         data{r}.I = Idist;
        
        results = [C;Idist];
        labels = [repmat([region(1),'c'],length(C),1); repmat([region(1),'i'],length(Idist),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        %plotvars = [plotvars; results,labels];
    
    else 
         %calc significance
        %[~,p] = ttest2(C,I);
        p = ranksum(C,I);
        pvals = [pvals;p];
        disp([region, ' pvalue = ', num2str(p)]);
        
        %plot
%         errorbar(r,mean(C),Cerr,'k.');
%         hold on
%         errorbar(r+0.25,mean(I),Ierr,'k.');
%         data{r}.C = C;
%         data{r}.I = I;
        results = [C;I];
        labels = [repmat([region(1),'c'],length(C),1); repmat([region(1),'i'],length(I),1)];
        results_all = [results_all;results];
        labels_all = [labels_all;labels];
        
    end
    
    
    %rrPowerData.(region) = C;
    %save([topDir,'AnalysesAcrossAnimals\rrPowerData'],'rrPowerData');
end

xticks([1.125 2.125 3.125 4.125]);
xticklabels(regions);

axis([0.75 3.75 -0.2 1.2]);
text(1.5, -.2, ['p = ', num2str(round(pvals(1),2))]);
text(2.5, -.2, ['p = ', num2str(round(pvals(2),2))]);
text(3.5, -.2, ['p = ', num2str(round(pvals(3),2))]);
%text(4.5, -.2, ['p = ', num2str(round(pvals(4),2))]);

ylabel('Respiratory Rhythm (7-8 Hz) Power');
figfile = [figDir, 'EEG\RespPower_performance'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile); 

figure
hold on
boxplot(results_all,labels_all,'Colors','br','symbol','k.');
xticks([1.5 3.5 5.5 7.5]);
xticklabels(regions);
text(1.5, 4, ['p = ', num2str(round(pvals(1),2))]);
text(3.5, 4, ['p = ', num2str(round(pvals(2),2))]);
text(5.5, 4, ['p = ', num2str(round(pvals(3),2))]);
%text(7.5, 7, ['p = ', num2str(round(pvals(4),2))]);
%axis([0 r+1.25 -0.1 1.1])
%text(1,0.16,['p = ', num2str(p)]);
%plot([0 r+1.25],[0 0],'k--');
ylabel('Respiratory Rhythm (7-8 Hz) Power');
%legend({'Correct', 'Incorrect'},'Location','East')

figfile = [figDir, 'EEG\RespPower_performance_boxplot'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);
