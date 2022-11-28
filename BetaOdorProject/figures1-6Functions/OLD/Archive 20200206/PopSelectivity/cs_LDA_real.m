%do the LDA separately for all animals and days, end up with a distribution
%for each day? use means?
%might not need to do PCA in this case, since number of selective cells on
%any given day should be <10
close all
clear
topDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
totalwin = [0 1];
regions = {'CA1','PFC'};
%win(1) = 0 - win(1); %make negative
binsize = 0.05;
bins = totalwin(1)+binsize:binsize:totalwin(2);
createvectors =1;
for r = 1:length(regions)
    region = regions{r};
    
    LDAvectors = [];
    mndist_all = [];
    
    if createvectors ==1
        
        for t = 1:length(bins)
            win = [0 bins(t)];
            
            mndist = [];
            daynum = 0;
            for a = 1:length(animals)
                animal = animals{a};
                animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
                runeps = cs_getRunEpochs(animDir, animal, 'odorplace');
                days = unique(runeps(:,1));
                
                for day = days'
                    disp(['Doing 0-',num2str(win(2)),'s ',animal,' day ', num2str(day)])
                    [leftcol,rightcol] = cs_columnVectors_day(animal,day,region,win,1);
                    col = [leftcol;rightcol];
                    lefttrials = size(leftcol,1);
                    righttrials = size(rightcol,1);
                    if size(col,2) < 3 %only use days with minimum number of selective cells
                        continue
                    
                    end
                    daynum = daynum+1;
                    Target = [zeros(lefttrials,1);ones(righttrials,1)];
                    LDAvectors{t}{daynum} = [col,Target];
                    %LDA column vector cell array:
                    %LDAvectors{bin size}{day number}(trials x cells)
                end
            end
        end
        save([topDir, 'AnalysesAcrossAnimals\LDAvectors_',region],'LDAvectors');
    else
        load([topDir, 'AnalysesAcrossAnimals\LDAvectors_',region]);
    end
    
    for t = 1:length(bins)
        win = [0 bins(t)];
        binvectors = LDAvectors{t};
        for d = 1:length(binvectors)
            col = binvectors{d}(:,1:end-1);
            Target = binvectors{d}(:,end);
            
            
            W = LDA(col,Target); %get the LD
            
            L = [ones(size(col,1),1) col] * W'; %project the data onto the LD
            
            mn= mean(mean(pdist2(L(1:lefttrials,:), L(lefttrials+1:end,:))));
            mndist = [mndist;mn];
            mndist_all{t} = mndist;
            %dists_cheby = max(abs(bsxfun(@minus, L(1:43,:), L(44:end,:))),[],2);
            
            %         figure
            %         plot(L(1:lefttrials,1),L(1:lefttrials,2),'r.')
            %         hold on
            %         plot(L(lefttrials+1:end,1),L(lefttrials+1:end,2),'b.')
            
            
            
            %         [signals,PC,V] = pca1(L');
            %         figure
            %         plot(signals(1,1:lefttrials),0,'r.')
            %         hold on
            %         plot(signals(1,lefttrials+1:end),0,'b.')
            %
            %         g1 = [repmat(1,lefttrials,1);repmat(2,righttrials,1)];
            %          p = anovan(signals(1,:)',g1)
            %
            %         L1 = [L(1:lefttrials,1);L(lefttrials+1:end,1)];
            %         L2 = [L(1:lefttrials,2);L(lefttrials+1:end,2)];
            %         %g1 = [repmat(1,size(L,1),1);repmat(2,size(L,1),1)];
            %         g1 = [repmat(1,lefttrials,1);repmat(2,righttrials,1)];
            %         %L1 = [L(1:lefttrials,1),L(1:lefttrials,2);L(lefttrials+1:end,1),L(lefttrials+1:end,2)];
            %
            %         p = anovan(L1,g1)
            %         p = anovan(L2,g1)
            
            %         p = anovan([L1;L2],{[g1;g1]})
            %
            %         [d,p,stats] = manova1(L,g1)
            
            
            
        end
        
    end
end
