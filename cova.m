function y = cova(x)
% coefficient of variation
  x = x(:);
  y = std(x) / mean(x);
end
