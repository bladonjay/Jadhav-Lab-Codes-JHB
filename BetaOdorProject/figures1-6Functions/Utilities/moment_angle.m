function m = moment_angle(a,order,dim)
%MOMENT_ANGLE Uncentered trigonometric moments of all orders.
%
%   M = MOMENT_ANGLE(A,ORDER) returns the (complex-valued) ORDER-th uncentered
%   trigonometric moment of the values in A. For vector input, M is a scalar.
%   For matrix input, M is a row vector containing the trigonometric moment of
%   each column of A. For N-D arrays, MOMENT_ANGLE operates along the first
%   non-singleton dimension.
%
%   MOMENT_ANGLE(A,ORDER,DIM) takes the moment along dimension DIM of A.
%
%   MOMENT_ANGLE returns uncentered moments. To get central moments, subtract
%   the sample mean. See also: MEAN_ANGLE, MINUS_ANGLE.
%
%Written by SMK, 2009 December 9
%

if ~isscalar(order) || ~isnumeric(order) || ~(order > 0) || ...
    ~(round(order) == order)
  error('ORDER must be a positive integer');
end

if (nargin < 3) || isempty(dim)
  % The output size for [] is a special case, handle it here.
  if isequal(a,[])
    m = nan(class(a));
    return; 
  end;
  % Figure out which dimension mean will work along
  dim = find(size(a) ~= 1,1,'first');
  if isempty(dim)
    dim = 1;
  end
end

m = mean(complex(cos(order*a),sin(order*a)),dim);

