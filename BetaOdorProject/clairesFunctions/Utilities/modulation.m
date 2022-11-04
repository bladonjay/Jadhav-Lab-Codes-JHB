function [moddepth peakphase] = modulation(phase, varargin) 
% function [moddepth peakphase] = modulation(phase, varargin) 
%	Computes the depth of modulation and the peak phase of the events whose
%	phase is given in the input variables.  Each phase is replaced with a
%	gaussian of area 1 and stdev pi/12, a histogram of phases is
%	constructed and the modulation depth is measured as peak / 
%       (peak + trough).
%
%	Note that phase can vary between either -pi and pi or 0 and 2pi
%
%       Options are 
%           'stdev', #  specifies the standard deviation of the gaussian
%           		function


stdev = pi/12;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'stdev'
                stdev = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    elses
        error('Options must be strings, followed by the variable');
    end
end

% if the phases are between 0 and 2 pi, shift them to be between -pi and pi
if (max(phase) > pi)
    negphase = find(phase > pi);
    posphase = find(phase <= pi);
    phase = [(phase(negphase) - 2*pi) ; phase(posphase)];
end

% histogram the phases with a very small bin
bins = -pi:stdev/100:pi;
count = histc(phase, bins);

% replication the counts 2 pi further up so that the smoothing works
bins = [bins (bins + 2*pi)];
count = [count ; count];

% smooth the count. Note that because we used stdev/100 as the bin increment,
% the gaussian is always 100 bins wide
g = gaussian(100, 600);

ps = smoothvect(count, g);

% compute the depth of modulation between 0 and 2 pi.  We ignore the edges
% because the smoothing is not circular

minbin = min(find(bins >= 0));
maxbin = max(find(bins <= 2 * pi));
bins = bins(minbin:maxbin);
ps = ps(minbin:maxbin);

% find the peak and the trough and calculate the depth of modulation
[mx i] = max(ps);
peakphase = bins(i);

[mn i] = min(ps);
troughphase = bins(i);
moddepth = mx / (mx + mn);
