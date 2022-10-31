[topDir,figDir] = cs_setPaths();
dataDir = [topDir, 'AnalysesAcrossAnimals\'];
regions = {'TC'};
 savetag = '';

for r = 1:length(regions)
    region = regions{r};
    

load([dataDir, 'PSDdata_', savetag,region])

fpass = PSDdata.fpass;
data = PSDdata.data;
psd = mean(data,2);
stepsize = (fpass(2)-fpass(1))/length(psd);
freqs = fpass(1):stepsize:fpass(2)-stepsize;

figure
plot(freqs,psd);
ylabel('Power');
xlabel('Frequency');

figfile = [figDir,'EEG\PSD_',savetag,region,'_',num2str(fpass(1)),'-',num2str(fpass(2)),'Hz'];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);
end