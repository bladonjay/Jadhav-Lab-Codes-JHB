%CS 4/4/2018
%edited 12/12/2019
%calculates PCA for left and right trials. Can either consider all trials
%individually, or just take mean firing rate of each cell across trials. 
%plots new vectors for PCs 1-3. Then performs a shuffle to determine at
%which bins the right and left vectors significantly diverge. 

function cs_PCA(region, win, binsize, selectiveonly)
%close all
%cs_PCA('CA1',[0 1], 0.05, 1)

%% --- Params 
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
[dataDir, figDir] = cs_setPaths();
figdir = [figDir,'PCA\']; 

digits(16);

win1 = 0 - win(1); %make negative

winstr = [(num2str(win1*1000)),'-',num2str(win(2)*1000),'ms'];

if selectiveonly == 1
    selstr = '_selectivecells';
else
    selstr = '';
end

%% --- Calculate
disp('creating column vectors')
[leftfr, rightfr, lefttrials, righttrials, cellinds] = cs_columnVectors(animals, region, win, binsize, selectiveonly, 0);

%allbins = [-0.5:binsize:1.5-binsize];
goodbins = [win1:binsize:win(2)-binsize]; 
%binind = lookup(goodbins,allbins);
%numtimebins = length(goodbins);

numtimebins = length(goodbins); 

    LRcombined = [leftfr, rightfr];
    
disp('calculating pca')   
    [coeff,score] = pca(LRcombined');
    
    new = score(:,1:2)*coeff(1:2,:);
    newL = new(1:numtimebins,:);
    newR = new((numtimebins+1):end,:);
    
    newL = smoothdata(newL,1);
    newR = smoothdata(newR,1);
    
    
%     
%     LRcombined_long = [leftfr, rightfr];
%     [coeff2,score2,latent] = pca(LRcombined_long');
    
    %new2 = score2(:,1:2)*coeff(1:2,:);
%     
%         new2 = score2(:,1:2)*latent;
% 
%     newL2 = new2(1:length(allbins),:);
%     newR2 = new2((length(allbins)+1):end,:);
%     
%     newL2 = smoothdata(newL2,1);
%     newR2 = smoothdata(newR2,1);
    
    
%% --- Plot

%% --- Calculate the difference

diffPC1 = abs(newL(:,1)' - newR(:,1)'); 
diffPC2 = abs(newL(:,2)' - newR(:,2)');

figure,
set(gcf, 'Position', [50,100,1500,600]);
subplot(1,2,1)
hold on
% plot(goodbins, newL(:,1),'r-','LineWidth',3)
% plot(goodbins, newR(:,1),'b-','LineWidth',3)
plot(goodbins,diffPC1,'b-','LineWidth',3)
vline(0,'k--')
title('PC1')
ylabel('Left-Right PC difference')
xlabel('Time (s)')

subplot(1,2,2)
hold on
% plot(goodbins, newL(:,2),'r-','LineWidth',3)
% plot(goodbins, newR(:,2),'b-','LineWidth',3)
plot(goodbins,diffPC2,'b-','LineWidth',3)
vline(0,'k--')
title('PC2')
ylabel('Left-Right PC difference')
xlabel('Time (s)')



%% --- Save
%filename = ['pcaResults_',region,'_',winstr,'.mat' ];
%save([dataDir,'AnalysesAcrossAnimals\',filename],'prob_bin','coeff','score');

figtitle1 = ['PCA_',region,'_',winstr,selstr];

figfile = [figdir,figtitle1];

saveas(gcf,figfile,'fig');
print('-djpeg', figfile);
print('-dpdf', figfile);

%keyboard;


%% ---- Shuffle -----%%

iterations = 1000;
diffPC1_shuff = [];
diffPC2_shuff = [];

for i = iterations:-1:1
    if rem(i,50) == 0
        iterationnum = 1000-i+1;
        disp(['Doing iteration ',num2str(iterationnum)]);
    end
    
    [shuff_leftfr, shuff_rightfr] = cs_shufflespiking_v2(lefttrials, righttrials, cellinds, binsize);
    
    shuff_LRcombined = [shuff_leftfr, shuff_rightfr];
    
    warning('off','all')
    [coeff,score,~,~,~] = pca(shuff_LRcombined');
    
    new = score(:,1:2)*coeff(1:2,:);
    newL = new(1:numtimebins,:);
    newR = new((numtimebins+1):end,:);
    
    newL = smoothdata(newL,1);
    newR = smoothdata(newR,1);
    
    diffPC1_shuff(i,:) = abs(newL(:,1) - newR(:,1))';
    diffPC2_shuff(i,:) = abs(newL(:,2) - newR(:,2))';
end



% prob_bin = [];
% for p = 1:length(diffPC1)
%     binPC1 = diffPC1_shuff(:,p);
%     prob_bin(1,p) = sum(diffPC1(p)>binPC1)/length(binPC1);
% 
%     binPC2 = diffPC2_shuff(:,p);
%     prob_bin(2,p) = sum(diffPC2(p)>binPC2)/length(binPC2);
% end

shuffmean1 = mean(diffPC1_shuff);
shuffmean2 = mean(diffPC2_shuff);

shuffstd1 = std(diffPC1_shuff);
shuffstd2 = std(diffPC2_shuff);


%% --- Plot significance

subplot(1,2,1)
%plot(goodbins, shuffmean1,'g', 'LineWidth',3)
patch([goodbins, fliplr(goodbins)], [shuffmean1+(2*shuffstd1), shuffmean1-(2*shuffstd1)],'k','FaceAlpha',0.2);


subplot(1,2,2)
%plot(goodbins, shuffmean2,'g', 'LineWidth',3)
patch([goodbins, fliplr(goodbins)], [shuffmean2+(2*shuffstd2), shuffmean2-(2*shuffstd2)],'k','FaceAlpha',0.2);




% plot_prob_bin = prob_bin;
% plot_prob_bin(prob_bin<0.99) = nan;
% plot_prob_bin(prob_bin>=0.99) = 0;
% 
% 
% subplot(1,2,1)
% plot(goodbins, plot_prob_bin(1,:),'g', 'LineWidth',3)
% 
% subplot(1,2,2)
% plot(goodbins, plot_prob_bin(2,:),'g', 'LineWidth',3)

%% --- Save
filename = ['pcaResults_',region,'_',winstr,'.mat' ];
save([dataDir,'AnalysesAcrossAnimals\',filename],'coeff','score');

figtitle1 = ['PCA_',region,'_',winstr,selstr];

figfile = [figdir,figtitle1];

saveas(gcf,figfile,'fig');
print('-djpeg', figfile);
print('-dpdf', figfile);


