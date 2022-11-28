clear
dataDir = cs_setPaths();
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%win = [0 0.1];
%win(1) = 0 - win(1); %make negative
%binsize = 0.05;
region = 'CA1';

nsamp = 50;
loaddata = 0;

wins = [0.1:0.1:1];

% if loaddata == 1
%     load([dataDir,'AnalysesAcrossAnimals\columnVectors_randsamp_CA1_0-1000ms.mat']);
%     col = [leftfr;rightfr];
% else 
vals = [];
for w = 1:length(wins)
    %win = [0 wins(w)];
    win = [1 0];
    %win = [0 1];
    
    col = cs_columnVectors_randsamp(animals,region,win,nsamp,1,1);
% end
%[coeff,score,~,~,explained] = pca(col);


%use PCs that explain 80% of data
%good = find(cumsum(explained)<=80);
% coeff = coeff(good,good);
% score = score(:,good);


%pcs = score(:,1:2)*coeff(1:2,:);

[signals,PC,V] = pca1(col');
signals = signals';
explained = 100*V/sum(V);
good = find(cumsum(explained)<=50);
signals = signals(:,good);

Target = [zeros(nsamp,1);ones(nsamp,1)];
W = LDA(signals,Target); %get the LD



% what do i use for ANOVA? how do I project the data onto the LD? 
L = [ones(100,1) signals] * W'; %this is the projection! columns = dimensions, Nrows should = number of observations

% W = LDA(col,Target);
% L = [ones(100,1) col] * W';

figure
        plot(L(1:50,1),L(1:50,2),'r.')
        hold on
        plot(L(51:100,1),L(51:100,2),'b.')
        
        
%reshape L
L1 = [L(1:50,1),L(51:100,1);L(1:50,2),L(51:100,2)];

%[p] = anova2(L1,50)

%the two dimensions are repeated, so really the result is one-dimensional.
%Can just use the first dimension for calculating ANOVA. 
p = anova1(L1(1:50,:))

%p = anova1(L1(51:100,:))
% figure
%         plot(L1(1:50,1),0,'r.')
%         hold on
%         plot(L1(51:100,1),0,'b.')
%         
%         plot(mean(L1(1:50,1)),0,'ro');
%         plot(mean(L1(51:100,1)),0,'bo');

%pval = stats.pval;
vals = [vals; p];

end

plot(wins,vals);

goodbins = [win(1):binsize:win(2)-binsize]; 
numtimebins = length(goodbins); 

new = score(:,1:2)*coeff(1:2,:);
    newL = new(1:numtimebins,:);
    newR = new((numtimebins+1):end,:);
    
    Input = [newL(1,:),newR(1,:)];
    %Target = [zeros(size(newL,2),1); ones(size(newR,2),1)];
    
    W = LDA(Input,[1;2]);