%1/20/20
%Plot nice PSTH and rasters for specified odor selective cells. 
%pulls data from cs_odorSelectivity_v2 results
%Plots psth on correct and incorrec trials seprately for visualization


%Good cells for paper:
% 3_3_30_1- - best
% 7_4_2_9 
% 5_1_4_1 
%PFC:
% 8_1_15_2 
% 7_5_9_1 
% 2_2_22_2 -- best


%PF cells
% PFC:
% 8_2_23_7
%% Params
clear
[topDir, figDir] = cs_setPaths;
close all
regions = {'CA1','PFC'};

savefig = 1;


%bins = -0.45:0.05:1;
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};

leftcolor = rgb('BrightRoyalBlue');
rightcolor = rgb('Crimson');

%% Get data and plot
for r = 1:length(regions)
    region = regions{r};
    
    load([topDir, 'AnalysesAcrossAnimals\selectivityData_',region]);
    bins = selectivityData.win;
    %upsample data, then smooth
      newbins = bins(1)+0.05:0.05/5:bins(end);
    left_correct = smoothdata(interp1(bins(1:end),selectivityData.psth_left_correct',newbins),'gaussian',15)';
    right_correct = smoothdata(interp1(bins(1:end),selectivityData.psth_right_correct',newbins),'gaussian',15)';
    left_incorrect = smoothdata(interp1(bins(1:end),selectivityData.psth_left_incorrect',newbins),'gaussian',15)';
    right_incorrect = smoothdata(interp1(bins(1:end),selectivityData.psth_right_incorrect',newbins),'gaussian',15)';
    
   
    
    for c = 1:size(left_correct,1)
        index = selectivityData.cellinds(c,:);
        animal = animals{index(1)};
        figure
        subplot(1,2,1), set(gcf,'Position',[500 300 1200 400]);
        plot(newbins,left_correct(c,:),'LineWidth',3, 'Color', leftcolor);
        hold on
        plot(newbins,right_correct(c,:),'LineWidth',3,'Color',rightcolor);
        plot([0 0], [0 10],'k:')
        xlabel('Time from Odor Onset (s)');
        ylabel('Firing Rate (Hz)');
       
        
        subplot(1,2,2)
        plot(newbins,left_incorrect(c,:),'--','LineWidth',3, 'Color', leftcolor);
        hold on
        plot(newbins,right_incorrect(c,:),'--','LineWidth',3,'Color',rightcolor);
        plot([0 0], [0 10],'k:')
        
        xlabel('Time from Odor Onset (s)');
        ylabel('Firing Rate (Hz)');
        
        figtitle = ['PSTH_',region,'_',animal,'_',num2str(index(2)),'_',num2str(index(3)),'_',num2str(index(4))];
        figfile = [figDir,'OdorSelectivity\',figtitle];
        
        print('-dpdf', figfile);
        print('-djpeg', figfile);
        
        
    end
    close all
end