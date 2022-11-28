function [timebase,rawcorr,corr_sm] = spiketrainxcorr(spikes1,spikes2,bin,tmax,sw1)
%               
% Inputs:       spikes1 is a n x 1 vector containing n time stamps of spikes 
%               of cell 1
%
%               spikes2 is a m x 1 vector containing m time stamps of spikes 
%               of cell 2
%
%               bin is the bin size to compute the cross-correlation (s) 15
%  
%               tmax is the range of the time lags, so the window is
%               [-tmax, tmax] 500
%
%               sw1 is the size of the smoothing window (s)
%               
%               
% Outputs:      timebase is the time (lag) axis, = -tmax:bin:tmax
%               
%               rawcorr is the raw cross-correlation
%               
%               corr_sm is the smoothed cross-correlation
%               
% Description:  Performs cross-correlation on two input spike trains
%               
% Reference:    Euston, D. R., Tatsuno, M., & McNaughton, B. L. (2007). Fast-forward 
%               playback of recent memory sequences in prefrontal cortex during sleep. 
%               science, 318(5853), 1147-1150.
%               
% Author:       Wenbo TANG
%               
% Date:         March 7, 2018

timebase = -tmax:bin:tmax;

% create gaussian smooth window
nstd=round(sw1/(timebase(2) - timebase(1)));
x = 1:nstd+1;
mu = mean(x);
g0 = exp( (-(x-mu).^2)./(2*(nstd^2)));
s = sum(g0);
g1 = g0/s;

% calculate raw cross-correlation
rawcorr = zeros(1,length(timebase));
for tt = 1:length(timebase)
    for j = 1:length(spikes1)
         temp =find(abs(spikes2 - spikes1(j)-timebase(tt)) < bin/2);
         rawcorr(tt) = rawcorr(tt) + sum(length(temp));
     end
end

% normalization
rawcorr = rawcorr./length(spikes2);

% smooth
sm = conv(rawcorr, g1);
filterlen = length(g1);
% the convolution returns a vector whose length is the sum of
% the length of the two arguments - 1
startp = round(filterlen / 2);
endp = startp + length(rawcorr) - 1;
corr_sm = sm(startp:endp);
timebase = timebase.*-1;% flip, for canonical visualization
