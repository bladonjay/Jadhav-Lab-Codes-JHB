% g = GAUSSIAN2(stdev, size)
%     creates a 2D gaussian matrix SIZE x SIZE with the given stdev centered
%	  in the middle of g

function [g] = gaussian2(stdev, size)

g = zeros(size);

m = 1 + (size - 1) / 2;
mean = [m m];


for i = 1:size
	for j = 1:size
		x = [i j];
		g(i,j) = 1 / (stdev * sqrt(2 * pi)) * exp(-.5 * (mean - x) * ...
							(mean - x)' / stdev^2);
	end
end

% renormalize 
g = g ./ sum(sum(g));
