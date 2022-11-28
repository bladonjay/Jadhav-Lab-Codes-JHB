function cs_boxplot(x,g,varargin)
%creates boxplots of data. 
%input X should be data in a 1d array
%input g should be a 1d array the same size as x, indicating the groupings
%of the data in x

barwidth = 0.25;
linecolor = 'k';
plotdata = 1; %plot individual data points
jitter = 0.5*barwidth;

if nargin > 2
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'barwidth'
            barwidth = varargin{option+1};
        case 'linecolor'
            linecolor = varargin{option+1};
        case 'plotdata'
            plotdata = varargin{option+1};
        case 'jitter'
            jitter = 0.1;
    end
end
end

groups = unique(g);
figure, hold on
for n = 1:length(groups)
    group = groups(n);
    %find data corresponding to group
    dat = x(g == group);
    q = quantile(dat,[0.25 0.50 0.75]);
    
    %plot raw data
    if plotdata == 1
        scatter(repmat(group,length(dat),1), dat, 35, 'filled', 'jitter', 'on', 'jitterAmount', 0.15);
    end    
    
    %plot the box
    plot([group-barwidth group-barwidth group+barwidth group+barwidth, group - barwidth],[q(1) q(3) q(3) q(1) q(1)],linecolor);
    plot([group-barwidth group+barwidth], [q(2) q(2)],linecolor);
    
    %remove outliers
    o = isoutlier(dat);
    newdat = dat(~o);
    
    %find non-outlier max and min
    mx = max(newdat);
    mn = min(newdat);
    
    %plot whiskers
    plot([group group], [q(1) mn], [linecolor,'--']);
    plot([group group], [q(3) mx], [linecolor,'--']);
    
    %plot caps
    plot([group-(barwidth*0.5) group+(barwidth*0.5)], [mx mx], linecolor)
    plot([group-(barwidth*0.5) group+(barwidth*0.5)], [mn mn], linecolor)
    
    xticks(groups)

    
end