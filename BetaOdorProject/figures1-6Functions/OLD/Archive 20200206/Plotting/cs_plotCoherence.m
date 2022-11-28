function cs_plotCoherence(cohDataFile)
%close all
%'Coherence_CA1-OB.mat'
load(cohDataFile);

%----- Params ----- %

str= 'Tapers2-3';
if contains(cohDataFile,'fulltrial')
    str = [str,'_fulltrial'];
end
data = coherenceData.data;
animals = coherenceData.animals;
region = coherenceData.region;
trigtypes = coherenceData.trigtypes;
freqband = coherenceData.freqband;
fpass = coherenceData.fpass;

timewindow = coherenceData.data(1).output{1, 1}(1).timewindow  ;
%----- Calculate Averages -----%

% figure(100), hold on
% figure(200), hold on
for g = 1:length(trigtypes) %for each trigtype
    trigtype = trigtypes{g};
    
    cohtoget = ['cohgram_',trigtype];
    zcohtoget = [freqband,'coh_',trigtype];
    
    zcoherence = [];
    % ----- Average over animals -----%
    for c = 1:length(data) 
        
        cohgramdata = data(c).output{1,1};
        
        for i = 1:length(cohgramdata) %sum coherograms from all tetrode pairs
            
            cohgram(:,:,i) = cohgramdata(i).(cohtoget);
%             if ~exist('cohgram')
%                 cohgram = cohgramdata(i).(cohtoget);
%                 N = 1;
%             else
%                 cohgram = cohgram + cohgramdata(i).(cohtoget);
%                 N = N +1; 
%             end
            
        end
        
%         zcoherence = [zcoherence ; ([cohgramdata(i).(zcohtoget)])];
        
    end

cohgram = mean(cohgram,3); 


    
% std = 4;
%     s = gaussian2(std,(2*std));
%     cohgram = filter2(s,cohgram,'valid'); %'Valid' tag excludes edges that would be messed up by filter
    
    cohgram = interp2(cohgram,3);
    
    figure, 
    colormap(jet);
    %colormap(hot);
    imagesc(timewindow, 1:size(cohgram,1), cohgram(1:end,:))
    set(gca,'YDir','normal')
    colorbar
    hold on
    yticks(0:size(cohgram,1)/10:size(cohgram,1))
    yticklabels(0:10:fpass(2))
    plot([0 0],[1 size(cohgram,1)],'k--', 'LineWidth', 1.5);
    %yticklabels([10:10:fpass(2)])

    figtitle1 = ['Cohgram_',region,'_',str];
    title(sprintf('%s%d',figtitle1), 'Interpreter', 'none'); %allows underscores in fig titles
    

    figdir = 'D:\Figures\Cohgrams\';
    %figdir = 'E:\Figures\Cohgrams\';
    figfile = [figdir,figtitle1];
    
    print('-djpeg', figfile);
    print('-dpdf', figfile);
    saveas(gcf,figfile,'fig');
    
    clear cohgram

% ----- Plot z-scored Coherence in frequency band for trigtype ----- %
%     meanZCoherence = mean(zcoherence, 1);
%     semZCoherence = std(zcoherence,1)/ (sqrt(size(zcoherence,1)));
%     
%     figure(100),
%     plot(timewindow, meanZCoherence)
%     
%         switch trigtype
%             case 'allTriggers'
%                 plot(timewindow, meanZCoherence,'k')  
%             case 'correctTriggers'
%                 plot(timewindow, meanZCoherence,'g')
%                 patch([timewindow fliplr(timewindow)], [meanZCoherence+semZCoherence fliplr(meanZCoherence-semZCoherence)], 'g','FaceAlpha',.5 );
%             case 'incorrectTriggers'
%                 plot(timewindow, meanZCoherence,'r')
%                 patch([timewindow fliplr(timewindow)], [meanZCoherence+semZCoherence fliplr(meanZCoherence-semZCoherence)], 'r','FaceAlpha',.5 );
%         end
%     
%     
%     
% %----- Bar graph of mean coherence within bins before and after odor ----%
% 
%     
%       %get only -400ms - 400ms
%       %timebins = [-.4:.2:.4];
%       timebins = timewindow;
%       
%         %coh = vertcat(cohgramdata.(cohtoget));
%         cohnewbins = zeros(size(zcoherence,1),length(timebins)-1);
%         for b= 1:length(timebins)-1
%         
%         timewindow = round(timewindow,2);
%         bin = find(timewindow>=timebins(b) & timewindow<=timebins(b+1));
%         
%         cohinbin = zcoherence(:,bin);
%         meancohbin = mean(cohinbin,2);
%         
%         cohnewbins(:,b) = meancohbin;
%         end
%         
%       figure(200)
%       %x = [-9.5 -8.5 -7.5 -3.5 -2.5 -1.5 1.5, 2.5, 3.5 7.5, 8.5, 9.5]; 
%     x_a = timebins;
%     binsize = timebins(2) - timebins(1);
%     x_good = x_a(2:end)-(binsize/2);
%       meanCohBins = mean(cohnewbins,1)';
%       semCohBins = (std(cohnewbins,1)/sqrt(size(cohnewbins,1)))';
%       
%       switch trigtype
%             case 'allTriggers'
%                 %x_good = x(1:3:end);
%                 bar(x_good,meanCohBins, 'c'), hold on 
%                 for b = 1:length(x_good)
%                     alpha = 0.05/(length(x_good));
%                     [h_a, p_a] = ttest(cohnewbins(:,b),0,alpha);
%                     stats_a(1:2,b) = [h_a, p_a];
%                     if h_a ==1
%                         plot(x_good(b), meanCohBins(b)+0.1, 'k*', 'MarkerSize',5)
%                     end
%                 end
%                 errorbar(x_good, meanCohBins,semCohBins,'LineStyle','none', 'Color','k');
%                 
%                % stats_a = [ h_a, p_a ];
%             case 'correctTriggers'
% %                 x_good = x(2:3:end);
% %                 bar(x_good,meanCohBins, 'g', 'BarWidth', 0.15), hold on
% %                 [h_c, p_c] = ttest2(cohnewbins(:,1), cohnewbins(:,3));
% %                  stats_c = [h_c, p_c];
%             case 'incorrectTriggers'
% %                 x_good = x(3:3:end);
% %                 bar(x_good,meanCohBins, 'r', 'BarWidth', 0.15), hold on
% %                 [h_i, p_i] = ttest2(cohnewbins(:,1), cohnewbins(:,3));
% %                 stats_i = [h_i, p_i];
%       end       
%       
%       
%       
%         
% 
%         
% %     bins = [-.4:.2:.4];
% %     for b = 1:length(bins)
% %         cohinbin = mean(beta((beta(:,1)>=bin(1) & beta(:,1)<=bin(2)),2));
% %     zcoherence
%     
% end
% figure(100)
% figtitle2 = [freqband,'Coherence_',region];
% plot([0 0], [-0.25 1.0], 'k--', 'LineWidth', 2);
% title(sprintf('%s%d',figtitle2), 'Interpreter', 'none'); 
% xlabel('Time from nosepoke trigger')
% ylabel(['Average z-scored ',freqband, ' coherence']);
% 
% figfile = [figdir,figtitle2];
% 
% 
% print('-djpeg', figfile);
% print('-dpdf', figfile);
% saveas(gcf,figfile,'fig');
% 
% 
% 
% coherenceData.stats_all = stats_a ;
% % coherenceData.stats_correct = stats_c;
% % coherenceData.stats_incorrect = stats_i;
% 
% figure(200)
% plot([0 0], [0 .5], 'k--')
%         
% %set(gca,'XTickLabels', [-.2 0 .2])
% xlabel('Seconds from Odor Onset');
% ylabel('Average Z scored beta coherence');
% 
% figtitle3 =[freqband, 'Coherence', region,'_binned'];
% title(sprintf('%s%d',figtitle3), 'Interpreter', 'none'); 
% 
% figdir = 'D:\Figures\EEG\';
% %figdir = 'E:\Figures\EEG\';
% figfile = [figdir,figtitle3];
% print('-dpng', figfile);
% print('-djpeg', figfile);
% saveas(gcf,figfile,'fig');
% 
% filename =  ['D:\OdorPlaceAssociation\AnalysesAcrossAnimals\Coherence_',region,'.mat'];
% %filename = ['E:\AnalysesAcrossAnimals\Coherence_',region,'.mat'];
% save(filename,'coherenceData');

end

