function [cx lags] = wincorr(x,y,nlags)
% returns cross correlation of x,y where y is smaller
% lag=0 corresponds to the x(end),y(end) lined up

% TODO multiply in frequency domain

  [nx nchannels] = size(x);
  ny = size(y,1);

  % if this holds, then inner loop can be vectorized
  assert(nx >= 2*ny)

  % mean-center both signals
  x = bsxfun(@minus, x, mean(x,1));
  y = bsxfun(@minus, y, mean(y,1));

  cx = nan(1+nlags,nchannels);
  jj = 1:ny;
  for lag = 0:nlags
    x_ = x(nx-ny+jj-lag,:);
    c = sum(y .* x_, 1);
    cx(end-lag,:) = c;
  end

  cx = bsxfun(@rdivide, cx, cx(end,:)); % normalize
  lags = -nlags:0;
end
