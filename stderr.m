function e = stderr(x, dim)
% standard error
  if nargin < 2
    dim = find(size(x) > 1, 1); % first nonsingleton dimension
  end
  e = std(x,0,dim) / sqrt(size(x,dim));
end
