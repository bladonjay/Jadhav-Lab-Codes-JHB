%cs_betaLearning

%compare beta power and coherence between early and late learning. 

close all
clear
figure
hold on
[topDir,figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1-PFC','CA1-OB','PFC-OB'};

for r = 1:length(regions)
    Pre = [];
    Post = [];
    Fam = [];

    region = regions{r};
    disp(['Doing ', region])
    for a = 1:length(animals)
        animal = animals{a};
        disp(['Starting animal ',animal]);
        animDir = ([topDir, animal,'Expt\',animal, '_direct\']);
        dm1 = cs_getRunEpochs(animDir, animal, 'novelodor');
        dm2 = cs_getRunEpochs(animDir, animal, 'novelodor2');
        daymatrix = sortrows([dm1;dm2],1);
        tetinfo = loaddatastruct(animDir, animal,'tetinfo');
        
        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
            npWins = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
            
            epochs = daymatrix(daymatrix(:,1) == day,2);
            load([animDir, animal, 'coherence',region,daystr]);
            
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                time = coherence{day}{ep}.time;
                data = coherence{day}{ep}.Coh;
                freq = coherence{day}{ep}.freq;
                
                if isempty(data)
                    continue
                end
                %find beta freq rows
                betarows = freq >= 20 & freq <=30;
                data = data(betarows,:);
                data = mean(data,1); %mean beta coherence
                
                prelearnwins = npWins{day}{ep}(ismember(odorTriggers{day}{ep}.allTriggers,odorTriggers{day}{ep}.prelearn),:);
                postlearnwins = npWins{day}{ep}(ismember(odorTriggers{day}{ep}.allTriggers,odorTriggers{day}{ep}.postlearn),:);
                
                for t = 1:size(prelearnwins)
                    pre = mean(data(isExcluded(time,prelearnwins(t,:))));
                    Pre = [Pre;pre];
                end
                for t = 1:size(postlearnwins)
                    post = mean(data(isExcluded(time,postlearnwins(t,:))));
                    Post = [Post;post];
                end
                
            end
        end
        %% Familiar
                
        daymatrix = cs_getRunEpochs(animDir, animal, 'odorplace');

        days = unique(daymatrix(:,1));
        for day = days'
            daystr = getTwoDigitNumber(day);
            odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',day);
            npWins = loaddatastruct(animDir, animal, 'nosepokeWindow',day);
            
            epochs = daymatrix(daymatrix(:,1) == day,2);
            load([animDir, animal, 'coherence',region,daystr]);
            
            for ep = epochs'
                epstr = getTwoDigitNumber(ep);
                
                time = coherence{day}{ep}.time;
                data = coherence{day}{ep}.Coh;
                freq = coherence{day}{ep}.freq;
                
                if isempty(data)
                    continue
                end
                %find beta freq rows
                betarows = freq >=20 & freq <=30;
                data = data(betarows,:);
                data = mean(data,1); %mean beta coherence
                
                famwins = npWins{day}{ep};

                for t = 1:size(famwins)
                    fam = mean(data(isExcluded(time,famwins(t,:))));
                    Fam = [Fam;fam];
                end

            end
        end
    end
    Pre = Pre(~isnan(Pre));
    Preerr = stderr(Pre);
    
    Post = Post(~isnan(Post));
    Posterr = stderr(Post);
    
    Fam = Fam(~isnan(Fam));
    Famerr = stderr(Fam);
    
    
    p = ranksum(Pre,Post);
        disp([region, ' pvalue = ', num2str(p)]);
    
    figure
    errorbar(1, mean(Pre), Preerr,'k.')
    hold on
    errorbar(2, mean(Post), Posterr,'k.')
    errorbar(3, mean(Fam), Posterr,'k.')
    switch region
        case 'CA1-PFC'
            axis([0 4 0 0.25])
            text(1,0,['p = ',num2str(p)]); 
        case 'PFC-OB'
            axis([0 4 -.2 0.1])
            text(1,0,['p = ',num2str(p)]); 
        case 'CA1-OB'
            axis([0 4 0 0.3])
            text(1,0,['p = ',num2str(p)]); 
    end
    xticks([1 2 3])
    xticklabels({'pre-learn','post-learn','familiar'})
    ylabel('Z-scored Beta Coherence');
    
    figfile = [figDir,'EEG\betaCoherenceLearning_',region];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
    
end