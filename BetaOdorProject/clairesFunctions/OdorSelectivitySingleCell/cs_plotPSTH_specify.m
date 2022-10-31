clear
[topDir,figDir] = cs_setPaths;
animal = 'CS39';
region = 'PFC';
%day = 1;
cell = [7 9 2];
day = cell(1);
runtype = 'noodor';


animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
spikes = loaddatastruct(animDir, animal, 'spikes',day);
odorTriggers = loaddatastruct(animDir, animal,'odorTriggers',day);
leftcolor = rgb('BrightRoyalBlue');
rightcolor = rgb('Crimson');

runeps = cs_getRunEpochs(animDir, animal, runtype,day);
runeps = (runeps(:,2));

win = [0.45 1];
binsize = 0.05;

correctleftspikes = [];
correctrightspikes = [];
incorrectleftspikes = [];
incorrectrightspikes = [];
                
for ep = 1:length(runeps)
    epoch = runeps(ep);
    
    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
        epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
        
        [correct_left, correct_right, incorrect_left, incorrect_right] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
        %triginds = sort([correct_left; correct_right]);
        
        trigs = odorTriggers{day}{epoch}.allTriggers;
        %% Find Spikes on each trial
        for t = 1:length(trigs)
            trigwin = [trigs(t)-win(1), trigs(t)+win(2)];
            winspikes = epspikes(epspikes > trigwin(1) & epspikes <= trigwin(2));
            bins = (trigs(t)-win(1):binsize:trigs(t)+win(2));
            binspikecount = histcounts(winspikes,bins);
            
            %classify into trial type
            if ismember(t, correct_left)
                correctleftspikes = [correctleftspikes; binspikecount];
            elseif ismember(t, correct_right)
                correctrightspikes = [correctrightspikes; binspikecount];
            elseif ismember(t, incorrect_left)
                incorrectleftspikes = [incorrectleftspikes; binspikecount];
            elseif ismember(t, incorrect_right)
                incorrectrightspikes = [incorrectrightspikes; binspikecount];
            end
            
        end
        
    end
    
end

    bins = -win(1):binsize:win(2);

cleftfr = mean(correctleftspikes,1)./binsize;
crightfr = mean(correctrightspikes,1)./binsize;

cleftfr_err = stderr(correctleftspikes)./binsize;
crightfr_err = stderr(correctrightspikes)./binsize;

seldat = (nansum(cleftfr) - nansum(crightfr))./(nansum(cleftfr) + nansum(crightfr));

ileftfr = mean(incorrectleftspikes,1)./binsize;
if isempty(ileftfr)
    ileftfr = zeros(1,length(bins)-1);
end
irightfr = mean(incorrectrightspikes,1)./binsize;
if isempty(irightfr)
    irightfr = zeros(1,length(bins)-1);
end
%ileftfr_err = stderr(ileftfr)./binsize;
%irightfr_err = stderr(irightfr)./binsize;
    g = 15;
    newbins = bins(1)+binsize:binsize/3:bins(end);
    left_correct = smoothdata(interp1(bins(2:end),cleftfr,newbins),'gaussian',g)';
    right_correct = smoothdata(interp1(bins(2:end),crightfr,newbins),'gaussian',g)';
    %left_cerr = smoothdata(interp1(bins(2:end),cleftfr_err,newbins),'gaussian',g)';
    %right_cerr = smoothdata(interp1(bins(2:end),crightfr_err,newbins),'gaussian',g)';

    
    left_incorrect = smoothdata(interp1(bins(2:end),ileftfr,newbins),'gaussian',g)';
    right_incorrect = smoothdata(interp1(bins(2:end),irightfr,newbins),'gaussian',g)';
    %left_ierr = smoothdata(interp1(bins(2:end),ileftfr_err,newbins),'gaussian',g)';
    %right_ierr = smoothdata(interp1(bins(2:end),irightfr_err,newbins),'gaussian',g)';
    
    maxy = round(max(max([left_correct;right_correct;left_incorrect;right_incorrect]))+1);
figure
        subplot(1,2,1), set(gcf,'Position',[500 300 1200 400]);
        plot(newbins,left_correct,'LineWidth',3, 'Color', leftcolor);
        hold on
        plot(newbins,right_correct,'LineWidth',3,'Color',rightcolor);
        plot([0 0], [0 maxy],'k:')
        xlabel('Time from Odor Onset (s)');
        ylabel('Firing Rate (Hz)');
       
        %patch([newbins,fliplr(newbins)],[left_correct+left_cerr; flipud(left_correct-left_cerr)]','k','FaceAlpha',0.5);
        %patch([newbins,fliplr(newbins)],[right_correct+right_cerr; flipud(right_correct-right_cerr)]','k','FaceAlpha',0.5);

        
        
        subplot(1,2,2)
        plot(newbins,left_incorrect,'--','LineWidth',3, 'Color', leftcolor);
        hold on
        plot(newbins,right_incorrect,'--','LineWidth',3,'Color',rightcolor);
        plot([0 0], [0 maxy],'k:')
        
        %patch([newbins,fliplr(newbins)],[left_incorrect+left_ierr; flipud(left_incorrect-left_ierr)]','k','FaceAlpha',0.5);
        %patch([newbins,fliplr(newbins)],[right_incorrect+right_ierr; flipud(right_incorrect-right_ierr)]','k','FaceAlpha',0.5);

        
        xlabel('Time from Odor Onset (s)');
        ylabel('Firing Rate (Hz)');
        
        figtitle = ['PSTH_',region,'_',animal,'_',num2str(cell(1)),'_',num2str(cell(2)),'_',num2str(cell(3))];
        figfile = [figDir,'OdorSelectivity\',figtitle];
        
        print('-dpdf', figfile);
        print('-djpeg', figfile);

%allSI = [allSI; cellinfo{cell(1)}{runeps(1)}{cell(2)}{cell(3)}.SI];

