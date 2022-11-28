function stats = rayleigh_test(a)
%RAYLEIGH_TEST Test isotropy of circular data against a unimodal alternative.
%
%   STATS = RAYLEIGH_TEST(A) tests the hypothesis that the sample of angles (in
%   radians) in the vector A are drawn from a uniform circular distribution,
%   versus the alternative that the sample is drawn from a unimodal
%   distribution. This test is not designed to infer multimodality! The return
%   struct STATS has the following fields:
%       n: number of elements in A
%   theta: sample estimate of mean
%       Z: test statistic (see section 4.3 of Fisher, 1993)
%       p: p-value
%
%   Reference:
%   [1] Fisher N.I. (1993) _Statistical Analysis of Circular Data_. Cambridge
%       University Press.
%
%Depends on:
%   MOMENT_ANGLE (written by SMK)
%
%Written by SMK, 2010 January 07.
%

if ~isvector(a) || ~isfloat(a) || ~isreal(a) || ...
    any(a == +Inf) || any(a == -Inf)
  error('A must be a vector of floating-point real non-infinite values');
end
nan_flag = isnan(a);
if any(nan_flag)
  warning('A contains NaN values. These will be treated as missing values');
  a(nan_flag) = [];
end
if any((a(:) < -pi) | (a(:) > pi))
  warning('some A values are outside the interval [-pi,+pi]');
end

m = moment_angle(a,1);
stats.n = numel(a);
stats.theta = angle(m);
stats.Z = stats.n * (m .* conj(m));
stats.p = exp(-stats.Z) * (1 + (2*(stats.Z) - (stats.Z).^2)/(4*stats.n) - ...
    (24*stats.Z - 132*(stats.Z).^2 + 76*(stats.Z).^3 - ...
    9*(stats.Z).^4)./(288*(stats.n).^2) );


