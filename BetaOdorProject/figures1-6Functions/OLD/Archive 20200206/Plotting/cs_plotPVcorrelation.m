function cs_plotPVcorrelation(region, lowPDIbin, highPDIbin)

%cs_plotPVcorrelation('CA1', 3, 24)

[topDir, figDir] = cs_setPaths;

load([topDir,'AnalysesAcrossAnimals\PV.mat']);

PV = PV.(region);
left = PV.left;
right = PV.right;

left1 = left(:,lowPDIbin);
right1 = right(:,lowPDIbin);


mn = mean([left1,right1],2);

%diffs2 = left2 - right2;
[~,ind] = sort(mn);

% diffs1 = left1 - right1;
% [~,ind] = sort(diffs1);

left1 = left1(ind);
right1 = right1(ind);

[CCbin,p] = corrcoef(left1, right1);

figure, hold on
plot(left1,'b.','MarkerSize',15);
plot(right1,'r.','MarkerSize',15);

figtitle = 'ExamplePVCorrelation_high';
figfile = [figDir,'Spiking\',figtitle];
    
    print('-dpdf', figfile);
    print('-djpeg', figfile);

left2 = left(:,highPDIbin);
right2 = right(:,highPDIbin);

mn = mean([left2,right2],2);

%diffs2 = left2 - right2;
[~,ind] = sort(mn);

left2 = left2(ind);
%right2 = right2(ind);

figure, hold on
plot(left2,'b.','MarkerSize',15);
plot(right2,'r.','MarkerSize',15);

[CCbin,p] = corrcoef(left2, right2);

figtitle = 'ExamplePVCorrelation_low';
figfile = [figDir,'Spiking\',figtitle];
    
    print('-dpdf', figfile);
    print('-djpeg', figfile);
end
