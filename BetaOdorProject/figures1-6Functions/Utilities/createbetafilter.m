low = 15;
high = 30;

srate = 1500;
d = designeegfilt(srate,low,high);
kernel = d;
betafilter.kernel = kernel;
betafilter.samprate = srate;
betafilter.descript = [num2str(low),'-',num2str(high),' Hz beta filter'];
filterfile = ['D:\DataAnalysis\usrlocal\filtering\betafilter.mat'];
save(filterfile, 'betafilter');