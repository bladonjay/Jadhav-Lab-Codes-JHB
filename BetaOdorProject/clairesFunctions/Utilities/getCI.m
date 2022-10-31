function CI = getCI(x,p)
% x is a vector, matrix, or any numeric array of data. NaNs are ignored.
% p is the confidence level (ie, 95 for 95% CI)
% The output is 1x2 vector showing the [lower,upper] interval values.
CI = prctile(x,abs([0,100]-(100-p)/2));