%cs_odorTriggeredMUA

%This function plots multiunit activity across all animals during the
%nosepoke window

[topDir,figDir] = cs_setPaths();

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

buffer = 0.5; %ms before NP 
win = 1; %length of window to look at after NP

regions = {'CA1','PFC'};

for r = 1:length(regions)
    region = regions{r};
    mua = [];
    for a = 1:length(animals)
        animal = animals{a};
        
        animDir = [topDir,animal,'Expt\',animal,'_direct\'];
       
        runepochs = cs_getRunEpochs(animDir, animal, 'odorplace');
        
        days = unique(runepochs(:,1));
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            disp(['Doing ',animal,' day ',daystr]);
            %load spikes, cellinfo, and npwindow for each day
            nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow', day);
            spikes = loaddatastruct(animDir, animal, 'spikes', day);
            cellinfo = loaddatastruct(animDir, animal, 'cellinfo', day);
            
            epochind = runepochs(:,1) == day;
            epochs = runepochs(epochind,2);
            
            for ep = 1:length(epochs)
                epoch = epochs(ep);
                filt = ['isequal($area,''',region,''') & (isequal($tag, ''accepted'') | isequal($tag, ''mua''))'];
                cells = evaluatefilter(cellinfo{day}{epoch},filt);
                trigtimes = nosepokeWindow{day}{epoch}(:,1);
                npwindows = [trigtimes-buffer, trigtimes+win];
                
                for c = 1:size(cells,1)
                    cell = cells(c,:);
                    if ~isempty(spikes{day}{epoch}{cell(1)}{cell(2)}.data)
                    sp = spikes{day}{epoch}{cell(1)}{cell(2)}.data(:,1);
                    
                    
                    for t = 1:length(trigtimes)
                        trig = trigtimes(t);
                        goodspikes = sp(logical(isExcluded(sp, npwindows(t,:)))); %1's signify that spikes fall within time periods
                        goodspikes = goodspikes - trig;
                        mua = [mua;goodspikes];
                    end
                    end
                end
            end
        end
         
    end
     bins = -buffer:(buffer+win)/100:win;
     [counts,edges] = histcounts(mua,bins,'Normalization','pdf');
    
     newbins = edges(1):(edges(2)-edges(1))/3:edges(end);

     newcounts = smoothdata(interp1(bins(1:end-1),counts,newbins),'gaussian',20)';
figure
patch([newbins, fliplr(newbins)], [newcounts; zeros(length(newcounts),1)], 'k')
%plot(newbins,newcounts);
%axis([-buffer win 0 1])
hold on
plot([0 0], [0.5 0.8],'r:','LineWidth',3);
ylim([0.5, 0.8])

ylabel('Spiking Probability')
xlabel('Time from odor onset');

figfile = [figDir,'Spiking\odorTriggeredMUA_',region];
%saveas(gcf,figfile,'fig');
print('-dpdf', figfile);
print('-djpeg', figfile); 
close all
end