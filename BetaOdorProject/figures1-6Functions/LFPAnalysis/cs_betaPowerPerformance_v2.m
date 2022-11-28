%cs_betaPowerPerformance

%compare beta power between correct and incorrect trials
%average beta power across nosepoke window, get distribution for each trial
%type

%can also downsample, since there are way fewer incorrect trials
close all
clear
figure
hold on
[topDir,figDir] = cs_setPaths();
%animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};

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
        
        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
            npWins = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
            
            epochs = daymatrix(daymatrix(:,1) == day,2);
            load([animDir, animal, 'speclow',region,daystr]);
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                time = spec{day}{ep}.time;
                try
                data = spec{day}{ep}.Spec';
                catch
                    data = spec{day}{ep}.Coh';
                end
                freq = spec{day}{ep}.freq;
                
                %find beta freq rows
                betarows = freq >= 20 & freq <=30;
                data = data(betarows,:);
                data = mean(data,1); %mean beta power
                
                
                [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{ep});
                correctwins = npWins{day}{ep}([correct_left;correct_right],:);
                incorrectwins = npWins{day}{ep}([incorrect_left;incorrect_right],:);
                
                %get beta power for each trial
                for t = 1:size(correctwins)
                    correctbeta = mean(data(isExcluded(time,correctwins(t,:))));
                    C = [C;correctbeta];
                end
                for t = 1:size(incorrectwins)
                    incorrectbeta = mean(data(isExcluded(time,incorrectwins(t,:))));
                    I = [I;incorrectbeta];
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
        p = sum(Idist >= mean(C)) /length(Idist);
        disp([region, ' pvalue = ', num2str(p)]);
        text(r-0.1,0.2,['p = ', num2str(p)]);

        errorbar(r, mean(C), stderr(C),'b.')
        hold on
        errorbar(r+0.25, mean(Idist), 2*std(Idist),'r.')
        
        %plot(r, mean(C),'k.');
        %axis([0 r+1.25 0.1 0.75])
%         data{r}.C = mean(C);
%         data{r}.I = Idist;
    else 
         %calc significance
        %[~,p] = ttest2(C,I);
        p = ranksum(C,I);
        disp([region, ' pvalue = ', num2str(p)]);
        
        %plot
        errorbar(r,mean(C),Cerr,'k.');
        hold on
        errorbar(r+0.25,mean(I),Ierr,'k.');
%         data{r}.C = C;
%         data{r}.I = I;
    end
    
    
    betaPowerData.(region) = C;
    save([topDir,'AnalysesAcrossAnimals\betaPowerData'],'betaPowerData');
end

axis([0 r+1.25 -0.1 1.1])
%text(1,0.16,['p = ', num2str(p)]);
%plot([0 r+1.25],[0 0],'k--');
ylabel('Z-scored Beta Power');
xticks(1:length(regions))
xticklabels(regions);
legend({'Correct', 'Incorrect'},'Location','East')

figfile = [figDir, 'Behavior\betaPower_performance'];
            
    print('-djpeg', figfile);
    print('-dpdf', figfile);
