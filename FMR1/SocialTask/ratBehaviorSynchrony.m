% assay rat behavioral syncrhony




[myfile,mydir]=uigetfile(':/C','Find your matfile with all the behavioral data');

% loading that file should get you ratinfo into your workspace
load(fullfule(mydir,myfile));


%% how about whether they move together?
%{
basically do they alternate wells in close proximity?
the question i guess is what is the distribution of latencies between
animal arrival times, and is this greater than chance?

I think the way to do this is to take the xcorr of these arrival times at
say 1-30 second lags, and then see if there are any peaks... then run
this same analysis after circshifting one.

lets turn the timestamps into ones and zeros, something like accumarray
each poke in each well in 200 msec bins say, concatenate them, and any
times between different wells, those are transitions and those are ones,
all others are zeros, nontransitions

the other question is whether this actually swaps from time to time...
maybe one animal leads and so theres a significant bump behind him, and
then sometimes he follows, so the bump is ahead of him.  It may be a good
idea to do this within each session.  The best result would be that the
distribution is multimodal, and each animal has their favorite lag

one thing to do would be to remove all arrivals immediately following a
reward, because then they're locked to whenever the rats finish eating

%}




% lets just get the cc for each and every session

% input parameters to the correlation and plotter
binsize=1/4; maxlag=120; % in seconds
timeshifts=(-maxlag:binsize:maxlag);
plotIT=0;


wb=waitbar(0);
for i=1:length(ratinfo)
    
    % gather the last sample before and first after an arm transition
    sessStart=min([ratinfo(i).ratsamples{1}.start(1) ratinfo(i).ratsamples{2}.start(1)]);
    % lock rat 1 samples to sess start
    ratinfo(i).ratsamples{1}.start=ratinfo(i).ratsamples{1}.start-sessStart;
    ratinfo(i).ratsamples{1}.end=ratinfo(i).ratsamples{1}.end-sessStart;
    ratinfo(i).ratsamples{2}.start=ratinfo(i).ratsamples{2}.start-sessStart;
    ratinfo(i).ratsamples{2}.end=ratinfo(i).ratsamples{2}.end-sessStart;


    % first accum array each arm
    fulllength=ceil([max(table2array(ratinfo(i).ratsamples{1}(end,1:2))) ...
        max(table2array(ratinfo(i).ratsamples{2}(end,1:2)))]/binsize);
    transmat=zeros(max(fulllength),2);
    for ra=1:2
        mytrans=ratinfo(i).ratsamples{ra};
        % get all departure ts, multiply by binsize, add 1 so no zeros
        departures=accumarray(round(mytrans.end(diff(mytrans.thiswell)~=0)/binsize)+1,1);
        arrivals=accumarray(round(mytrans.start([0; diff(mytrans.thiswell)]~=0)/binsize)+1,1);
        
        for t=1:length(departures)-1
            if departures(t)==1 && length(arrivals)>t
                nextarrival=find(arrivals(t:end)==1,1,'first');
                if ~isempty(nextarrival)
                    transmat(t:t+nextarrival-1,ra)=1;
                end
            end
        end
    end
    
    if plotIT
        figure; plot((1:size(transmat))*binsize,transmat(:,1));
        hold on; plot((1:size(transmat))*binsize,transmat(:,2)+2); ylim([-1 4]);
    end
    myCmap=lines(4);
    % right direction means shift second input foreward
    samplecorr=xcorr(transmat(:,1),transmat(:,2),maxlag/binsize);
    ratinfo(i).samplecorr=samplecorr;
    %subplot(2,2,ra);
    

    % right is first vector to the right, e.g. rat 1 past, rat 2 future.
    % so if theres a dip in corr in rat 1 past and rat 2 future, that means
    % after rat 2 moves, rat 1 moves.
    bootcorr=[];
    % now for a bootstrap
    for bt=1:200
        myshift= max([randi(length(transmat)) maxlag*2/binsize]);
        myboot=circshift(transmat(:,1),randi(length(timeshifts)));
        bootcorr(bt,:)=xcorr(myboot,transmat(:,2),maxlag/binsize);
    end
    ratinfo(i).bootcorr=[prctile(bootcorr,1)' prctile(bootcorr,99)'];
    if plotIT==1
        figure;
        b=bar(timeshifts,samplecorr,1,'EdgeColor','none','FaceColor',myCmap(ra,:));
        hold on;
        b(2)=plot(timeshifts,prctile(bootcorr,1),'k--');
        plot(timeshifts,prctile(bootcorr,99),'k--');
    end
    waitbar(i/length(ratinfo),wb,'Running boots');
end
close(wb);



%%
% so i want to see if the rats behave the same way with the same partners,
% so this concatenates sessions for each genotype.

%       *** this does not include fx-wt pairs as of now***

% lets just gather one pair, say 201 and 204, the controls
cohort=[1 1 3 3];
mypair=[1 4;  2 3; 1 3; 2 4];
mytitle={'c1, ctrl-ctrl','c1 fx-fx','c3 ctrl-crtl','c3 fx-fx'};
%figure;
for i=1:4
subplot(2,2,i); 


thiscohort=ratinfo([ratinfo.cohortnum]==cohort(i));
mysess=thiscohort(cellfun(@(a) a(1)==mypair(i,1) && a(2)==mypair(i,2), {thiscohort.ratnums}));
flipsess=thiscohort(cellfun(@(a) a(1)==mypair(i,2) && a(2)==mypair(i,1), {thiscohort.ratnums}));


allsess=[cell2mat({mysess.samplecorr})'; fliplr(cell2mat({flipsess.samplecorr})')];
    
lownull=cell2mat(cellfun(@(a) a(:,1), {mysess.bootcorr},'UniformOutput',false))';
highnull=cell2mat(cellfun(@(a) a(:,2), {mysess.bootcorr},'UniformOutput',false))';
timeshifts=(-maxlag:binsize:maxlag);


mp=plot(timeshifts,nanmean(allsess));
% make my own boots here

allboot=[];
for bt=1:100
    bootshift=randi(length(timeshifts),size(allsess,1));
    bootrows=nan(size(allsess));

    for k=1:size(allsess,1)
        bootrows(k,:)=circshift(allsess(k,:),bootshift(k));
    end
    allboot(bt,:)=nanmean(bootrows);
end

hold on;
mp(2)=plot(timeshifts,prctile(allboot,99),'k--');
plot(timeshifts,prctile(allboot,1),'k--');
axis tight;
title(mytitle{i});
%xlabel(sprintf('Seconds offset \n Rat %d follows      Rat %d leads',mypair(1),mypair(1))); 
ylabel(sprintf('Correlation in movement \n behavior between rats'));
% it looks like their movement speed is anticorrelated...
% this is true for all the animals, so lets just show the fx-fx guys and
% the ctrl=ctrl pairs


%title('Cohort B, FX-FX pair');
xlabel(sprintf('Rat 1 follows                Rat 1 leads \n Seconds offset in Behavior')); 
legend(mp,'Real Data','Bootstrap 99% CI');
end

