signal = loadcached('data/ARN038_grid_block50_stim017.mat','signal');
[nch m] = size(signal);
npairs = nchoosek(nch,2);

% populate lookup table
lookup = nan(nch,nch);
lookup(1,2) = 0;
for i = 1:nch, for j = i+1:nch, lookup(i,j) = max(lookup(:))+1; end, end

% create time shifted matrices
X  = signal(:,1:end-1);
X_ = signal(:,2:end);

A = X_ * pinv(X);

for tt = 1:10
  norm(X(:,tt+1) - A * X(:,tt)) / norm(X(:,tt))
end

return

X = nan(npairs,n);
for i = 1:nch
  for j = 1:nch
  end
end
