function pval = cs_bootstrap(sample1, sample2, num)
numX = length(sample1);
numY = length(sample2);

diff_means = zeros(num,1);

pairs = [sample1, sample2];
for i = 1:num
    
    draw_X = datasample(sample1, numX);
    draw_Y = datasample(sample1, numY);
    
    diff_means(i) = mean(draw_X) - mean(draw_Y);
    mnx(i) = mean(draw_X);
    mny(i) = mean(draw_Y);
    
end

figure,
histogram(diff_means)

stdev = std(diff_means);

pval = mean(mnx > mny);