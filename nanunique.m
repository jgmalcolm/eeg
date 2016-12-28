function u = nanunique(x)
% omit redundant NaN output (only show once if present)
  isrow_ = isrow(x);
  x = x(:);
  ispresent = any(isnan(x));
  u = unique(x(~isnan(x)));
  if ispresent
    u(end+1) = nan;
  end
  if isrow_
    u = u';
  end
end
