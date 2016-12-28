function show_r2_rho(x,y,y_hat,unit)

  xx = unique(x);
  fprintf('\n%5s  %10s  %10s\n', unit, 'r2', 'rho')
  for i = 1:numel(xx)
    yy = y(x == xx(i));
    yy_hat = y_hat(x == xx(i));
    sse = sum((yy - yy_hat).^2);
    ss = sum((yy - mean(yy)).^2);
    r2 = 1-sse/ss;
    rho = corr(yy,yy_hat);
    fprintf('%5g  %10.5f  %10.5f\n', ...
            xx(i), r2, rho)
  end
end
