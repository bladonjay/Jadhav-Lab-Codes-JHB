function [m, err] = cs_sextiles(x,y,varargin)

if nargin > 2
    disp('test')
    n = varargin{1};
else
    n = 6;
end

[x,ind] = sort(x);
y = y(ind);

remainder = rem(length(x),n);
x_tmp = x(1:end-remainder);
y_tmp = y(1:end-remainder);
num = length(y_tmp)/n; %number of samples per bin

tmp = y_tmp;
sextiles = [];
for i = 1:n
    sextiles(:,i) = tmp(1:num);
    tmp(1:num) = [];
end

m = mean(sextiles);
err = std(sextiles)/sqrt(size(sextiles,1));
 errorbar(1:n,m,err)