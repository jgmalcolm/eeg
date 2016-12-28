function x = notnan(x)
  x(isnan(x)) = [];
end
