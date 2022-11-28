%CS 10/30/2018

%calculates PCA for left and right trials. Can either consider all trials
%individually, or just take mean firing rate of each cell across trials. 
%plots time x PC1 x PC2 for Right vs Left, calculates distance between 3d trajectories.
%Compares to shuffled data to calc significance. 

function cs_PCA_3d(region, win, binsize, selectiveonly)
%close all
%cs_PCA_3d('CA1',[0 1], 0.05, 1)

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

leftcolor = rgb('BrightRoyalBlue');
rightcolor = rgb('Crimson');

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
    [coeff,score,~,~,~] = pca(LRcombined');
    
    new = score(:,1:2)*coeff(1:2,:);
    newL = new(1:numtimebins,:);
    newR = new((numtimebins+1):end,:);
    
    newL = smoothdata(newL,1);
    newR = smoothdata(newR,1);
    
    realL = newL;
    realR = newR;
    
    
%% --- Plot
% figure,
% plot3(goodbins, newL(:,1), newL(:,2),'r-','LineWidth',3); 
% hold on
% plot3(goodbins, newR(:,1), newR(:,2),'b-','LineWidth',3);
% 
% grid on
% xlabel('Time from odor onset (seconds)');
% ylabel('PC1');
% zlabel('PC2');

%% --- Save
%filename = ['pcaResults_',region,'_',winstr,'.mat' ];
%save([dataDir,'AnalysesAcrossAnimals\',filename],'prob_bin','coeff','score');

% figtitle1 = ['PCA3d_',region,'_',winstr,selstr];
% 
% figfile = [figdir,figtitle1];
% 
% saveas(gcf,figfile,'fig');
% print('-djpeg', figfile);
% print('-dpdf', figfile);

%keyboard;

%% --- Calculate the difference

diffPC = sqrt(goodbins'.^2+(newL(:,1)-newR(:,1)).^2+(newL(:,2) - newR(:,2)).^2);
% figure(2),
% plot(goodbins,diffPC)

%% ---- Shuffle -----%%

iterations = 1000;
diffPC_shuff = [];

for i = iterations:-1:1
    if rem(i,20) == 0
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
    
    diffPC_shuff(i,:) = sqrt(goodbins'.^2+(newL(:,1)-newR(:,1)).^2+(newL(:,2) - newR(:,2)).^2);
end

prob_bin = [];
for p = 1:length(diffPC)
    binPC = diffPC_shuff(:,p);
    prob_bin(1,p) = sum(diffPC(p)>binPC)/length(binPC);

end

% shuffmean = mean(diffPC_shuff,1);
% 
% prc95 = prctile(diffPC_shuff,95,1);
% prc05 = prctile(diffPC_shuff,5,1);



%% --- Plot significance
% plot_prob_bin = prob_bin;
% plot_prob_bin(prob_bin<0.99) = nan;
% plot_prob_bin(prob_bin>=0.99) = 0;

nonsig = find(prob_bin< 0.99);
sig = find(prob_bin(nonsig(end)+1:end) >= 0.99)+nonsig(end);
sig = [nonsig(end),sig];

figure,
% plot3(  goodbins(nonsig), realL(nonsig,1),  realL(nonsig,2),'k-','LineWidth',4); 
% hold on
% plot3( goodbins(nonsig), realR(nonsig,1),   realR(nonsig,2),'k-','LineWidth',4);

plot3(  goodbins(nonsig), realL(nonsig,1),  realL(nonsig,2),'Color',leftcolor,'LineWidth',4); 
hold on
plot3( goodbins(nonsig), realR(nonsig,1),   realR(nonsig,2),'Color',rightcolor,'LineWidth',4);


plot3(   goodbins(sig), realL(sig,1), realL(sig,2),'Color',leftcolor,'LineWidth',4); 
plot3(   goodbins(sig),  realR(sig,1), realR(sig,2),'Color',rightcolor,'LineWidth',4);



ax = gca;
ax.LineWidth = 1.5;

view(-11.6, 12.72);
grid on
xlabel('Time from odor onset (seconds)');
ylabel('PC1');
zlabel('PC2');

% nonsigtime = goodbins(nonsig);
% sigtime = goodbins(sig);
% 
% nonsigdata = 



% plot3(goodbins, plot_prob_bin,plot_prob_bin,'g', 'LineWidth',3)

% figure(2) 
% hold on
% patch([goodbins, fliplr(goodbins)], [shuffmean-prc05,fliplr(shuffmean)+fliplr(prc95)], 'k','FaceAlpha',0.2,'EdgeColor','none')
%     p1 = plot(goodbins, shuffmean+prc95,'k-');
%     p1.Color(4) = 0.25;
%     p2 = plot(goodbins, shuffmean-prc05,'k-');
%     p2.Color(4) = 0.25;
%     plot(goodbins, shuffmean,'k-');


%% --- Save
%filename = ['pcaResults3d_',region,'_',winstr,'.mat' ];
%save([dataDir,'AnalysesAcrossAnimals\',filename],'prob_bin','coeff','score');

figtitle1 = ['PCA3dSignificance_',region,'_',winstr,selstr];

figfile = [figdir,figtitle1];

saveas(gcf,figfile,'fig');
print('-djpeg', figfile);
print('-dpdf', figfile);



