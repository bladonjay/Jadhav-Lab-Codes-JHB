%cs_cfc

%get cfc between RR and beta using method from Tort et al 2010
close all
clear
[topDir, figDir] = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC','OB'};
createdata = 1;
[A,B] = meshgrid(1:length(regions),1:length(regions));
pairInds=reshape(cat(2,A',B'),[],2);
regionpairs = regions(pairInds);

Nphasebins = 18; %20 degree bins
bins = -pi:(2*pi/Nphasebins):pi;
figure
for r = 1:size(regionpairs,1)
    r1 = regionpairs{r,1}; %RR
    r2 = regionpairs{r,2}; %beta
    
    if createdata == 1
        allphaseamp = [];
        alldecisiontimes = [];
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
                
                %get np off times after each high beta period
                
                wins = nosepokeWindows{day}{epoch};
                
                highbeta = highBeta{day}{epoch}.OB;
                %use only first 100ms of each beta window, controls for
                %more data on longer NPs
                highbeta(:,2) = highbeta(:,1)+0.1; 
                
                wins = wins(lookup(highbeta(:,2),wins(:,2)),:);
                wins = wins(:,2) - wins(:,1);
                
                %
                for t = 1:size(highbeta,1)
                    inds = isExcluded(lfptimes,highbeta(t,:));
                    
                    RRph = RRphase(inds);
                    bAmp = betaamp(inds);
                    allphaseamp{end+1} = [RRph,bAmp];
                    
                end
                alldecisiontimes = [alldecisiontimes;wins];
            end
        end
        
        
        save([topDir,'AnalysesAcrossAnimals/cfcData_',r1,'RR-',r2,'beta'],'allphaseamp','alldecisiontimes');
    else
        load([topDir,'AnalysesAcrossAnimals/cfcData_',r1,'RR-',r2,'beta']);
    end
    %get highest and lowest decision latencies
    [~,sortinds] = sort(alldecisiontimes);
    num = fix(0.25*length(sortinds));
    iterations = 1000;
    %% low
    MIdist_low = [];
    for i = 1:iterations
       
    lowlatency_time = alldecisiontimes(sortinds(1:num));
    lowlatency_phaseamp = allphaseamp(sortinds(1:num));
    
    samp = datasample(1:length(lowlatency_time), 500);
    
    lowlatency_time = lowlatency_time(samp);
    lowlatency_phaseamp = lowlatency_phaseamp(samp);
    
    %calculate mean phase-amplitude coupling
    phaseamp = cell2mat(lowlatency_phaseamp');
    binassign = discretize(phaseamp(:,1),bins);
    nans = isnan(binassign);
    phaseamp(nans,:) = [];
    binassign(nans) = [];
    
    meanPhAmp = accumarray(binassign,phaseamp(:,2),[],@mean);
    
    %normalize
    sumAmp = sum(meanPhAmp);
    normPhAmp = meanPhAmp./sumAmp;
    
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
    MIdist_low =[MIdist_low;MI];
    end
    %% high
    
    MIdist_high = [];
    for i = 1:iterations
    highlatency_time = alldecisiontimes(sortinds(end-num:end));
    highlatency_phaseamp = allphaseamp(sortinds(end-num:end));
    
        samp = datasample(1:length(lowlatency_time), 500);

        highlatency_time = highlatency_time(samp);
        highlatency_phaseamp = highlatency_phaseamp(samp);
        
    %calculate mean phase-amplitude coupling
    phaseamp = cell2mat(highlatency_phaseamp');
    binassign = discretize(phaseamp(:,1),bins);
    nans = isnan(binassign);
    phaseamp(nans,:) = [];
    binassign(nans) = [];
    
    meanPhAmp = accumarray(binassign,phaseamp(:,2),[],@mean);
    
    %normalize
    sumAmp = sum(meanPhAmp);
    normPhAmp = meanPhAmp./sumAmp;
    
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
    MIdist_high = [MIdist_high;MI];
    end
    subplot(3,3,r)
    errorbar([1 2],[mean(MIdist_low) mean(MIdist_high)],[2*std(MIdist_low) 2*std(MIdist_high)],'k.')
    xticks([1 2]);
    xticklabels({'low latency','high latency'})
    xlim([0 3])
    ylabel('MI');
    title([r1,'RR-',r2,'beta'])
    
    [~,lower] = min([max(MIdist_low),max(MIdist_high)]);
    %what % of lower distribution overlaps with higher distribution?
    p1 = sum(MIdist_low >= min(MIdist_high))/length(MIdist_low);
    p2 = sum(MIdist_high >= min(MIdist_low))/length(MIdist_high);
    p = min([p1 p2]);
    text(1, 0, ['p = ',num2str(round(p,2,'significant'))]); 
    
end

figfile = [figDir, 'EEG/cfc_decisionTime'];
print('-djpeg', figfile);
print('-dpdf', figfile);