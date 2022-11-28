%cs_cfc

%get cfc between RR and beta using method from Tort et al 2010
close all
clear
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};

[A,B] = meshgrid(1:length(regions),1:length(regions));
        pairInds=reshape(cat(2,A',B'),[],2);
        regionpairs = regions(pairInds);
        
        Nphasebins = 18; %20 degree bins
        bins = -pi:(2*pi/Nphasebins):pi;
        
for r = 1:length(regionpairs)
    r1 = regionpairs{r,1}; %RR
    r2 = regionpairs{r,2}; %beta
    
    allphaseamp = [];
    
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        dayepochs = cs_getRunEpochs(animDir,animal,'odorplace');
        
        nosepokeWindows = loaddatastruct(animDir, animal,'nosepokeWindow',dayepochs(:,1)');
        load([animDir,animal,'highBeta']);
        
        for d =1:size(dayepochs,1)
            day = dayepochs(d,1);
            epoch = dayepochs(d,2);
            
            %use tets with the most cells in each region
            tet1 = cs_getMostCellsTet(animal, day, epoch, r1);
            tet2 = cs_getMostCellsTet(animal, day, epoch, r2);
            
            RR = loadeegstruct(animDir,animal,'resp',day, epoch,tet1);
            RRphase = double(RR{day}{epoch}{tet1}.data(:,2))/10000;
            beta = loadeegstruct(animDir, animal, 'beta',day,epoch,tet2);
            betaamp = double(beta{day}{epoch}{tet2}.data(:,3));
            
            %zscore amplitude
            betaamp = (betaamp-mean(betaamp))/std(betaamp);
            
            lfptimes = geteegtimes(RR{day}{epoch}{tet1});
            %np = nosepokeWindows{day}{epoch};
            np = highBeta{day}{epoch}.OB;
            
            
            inds = isExcluded(lfptimes,np);
            
            RRphase = RRphase(inds);
            betaamp = betaamp(inds);
            allphaseamp = [allphaseamp;[RRphase,betaamp]];
            
        end
    end
    
    binassign = discretize(allphaseamp(:,1),bins);
    nans = isnan(binassign);
    allphaseamp(nans,:) = [];
    binassign(nans) = [];
    
    meanPhAmp = accumarray(binassign,allphaseamp(:,2),[],@mean);
    
    sumAmp = sum(meanPhAmp);
    normPhAmp = meanPhAmp./sumAmp;
    newbinedges = [bins - pi, bins(2:end) + pi];
    newbins = newbinedges(2:end)-((bins(2)-bins(1))/2);
    figure
    subplot(1,2,1)
    patch([newbins, fliplr(newbins)],[zeros(1,length(newbins)), repmat(normPhAmp',1,2)], 'k')
    %plot(newbins,repmat(normPhAmp',1,2))
    ylim([0 0.08])
    xlim([newbins(1) newbins(end)])
    xticks([newbins(1)  -pi 0 pi newbins(end)])
    xticklabels({'-2\pi','-\pi','0','\pi','2\pi'})

    xlabel('RR phase')
    ylabel('Normalized Beta Amplitude')
    
    
    %Kullback-Leibler distance
    Q = ones(length(normPhAmp),1);
    Q = Q ./sum(Q); %normalize so sums to 1. normPhAmp should already be normalized. 
    if sum(normPhAmp) ~= 1
        P = normPhAmp ./sum(normPhAmp);
    else
        P = normPhAmp;
    end
    tmp = P.*log(P./Q);
    dist = sum(tmp);
    
    %modulation index from Tort et al 2010
    MI = dist/log(length(P));
    
    
    %get statistical significance by shuffling ("surrogate control")
    MI_shuff= [];
    iterations = 500;
    for i = 1:iterations
        shuffamp = allphaseamp(randperm(length(allphaseamp)),2);
        shuffmn = accumarray(binassign,shuffamp,[],@mean);
        shuffsum = sum(shuffmn);
        shuffnorm = shuffmn./shuffsum;
        
        %Kullback-Leibler distance
    Q = ones(length(shuffnorm),1);
    Q = Q ./sum(Q); %normalize so sums to 1. normPhAmp should already be normalized. 
    if sum(shuffnorm) ~= 1
        P = shuffnorm ./sum(shuffnorm);
    else
        P = shuffnorm;
    end
    tmp = P.*log(P./Q);
    dist = sum(tmp);
    
    %modulation index from Tort et al 2010
    mi = dist/log(length(P));
    MI_shuff(i,1) = mi;
    end
    subplot(1,2,2)
    hold on
    bar([1 2], [mean(MI_shuff),MI]);
    numsig = 0.05*length(MI_shuff);
    sorted = sort(MI_shuff,1,'descend');
    thresh = sorted(numsig+1);
    plot([0 3], [thresh thresh],'k--')
    xticks([1 2])
    xticklabels({'Shuffled','Original'})
    ylabel('Modulation Index');
    
    p = sum(MI_shuff > MI)/length(MI_shuff);
    text(1, thresh, ['p = ',num2str(round(p,2,'significant'))]);
    figfile = [figDir,'EEG\CFC_',r1,'resp-',r2,'beta'];
    print('-djpeg',figfile);
    print('-dpdf',figfile);
    
end
    