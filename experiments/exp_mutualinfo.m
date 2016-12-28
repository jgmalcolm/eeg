clear
addpath('~/src/mutualinfo/')
addpath('~/src/Network-Controllability-Diagnostics/')
signal = loadcached('data/ARN038_grid_block50_stim017.mat','signal');
[nch m] = size(signal);
npairs = nchoosek(nch,2);

% populate lookup table
lookup = nan(nch,nch);
lookup(1,2) = 0;
for i = 1:nch, for j = i+1:nch, lookup(i,j) = max(lookup(:))+1; end, end

% center the signal
% signal = bsxfun(@minus, signal, mean(signal,2));
signal = zscore(signal,0,2);

A = zeros(nch); % diag==0 by definition
for i = 1:nch
  s1 = flat(signal(i,:));
  for j = i+1:nch
    s2 = flat(signal(j,:));
    [A(i,j) A(j,i)] = deal(mutualinfo(s1,s2));
    % [A(i,j) A(j,i)] = deal(corr(s1,s2,'type','Spear'));
  end
end
A

% [~,deg_rank] = sort(sum(A,1));
deg_rank = degTieredVals(A);
avg_ctrl = averMeas(A);
mod_ctrl = moduMeas(A);


chs = cell(8,1);
for i = 1:4, chs{i} = sprintf('CA3(%d)',i); end
for i = 5:8, chs{i} = sprintf('CA1(%d)',i); end

h = plot(1:8, avg_ctrl(deg_rank), 'r.--', ...
         1:8, mod_ctrl(deg_rank), 'b.--', ...
         'MarkerSize',30, 'LineWidth',3);
set(gca, 'FontSize', 16, ...
         'XTickLabel', chs(deg_rank), ...
         'LineWidth',1)
xlabel('Rank of mutual information')
ylabel('Controllability')
title('Rank of mutual information vs Average controllability')
legend(h, 'Average controllability', 'Modal controllability', 'Location','Best')
print('-dpng','-r200','results/rank_avg_controllability.png')

return

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
