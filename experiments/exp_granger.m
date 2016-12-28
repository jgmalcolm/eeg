clear
addpath('~/src/Network-Controllability-Diagnostics/')
addpath('~/src/mvgc_v1.0/')
fs = 6103.52;
signal = loadcached('data/ARN038_block50');
[nch m] = size(signal);

% center the signal
% signal = bsxfun(@minus, signal, mean(signal,2));
% signal = zscore(signal,0,2);

momax = 25;
acmaxlags = 1000;

[aic bic moaic mobic] = tsdata_to_infocrit(signal, momax, 'LWR');
[A SIG] = tsdata_to_var(signal, moaic);
assert(~isbad(A), 'VAR estimation failed')
[G info] = var_to_autocov(A, SIG, acmaxlags);
var_info(info,true)
F = autocov_to_pwcgc(G);
assert(~isbad(F,false), 'GC calculation failed')
F(1:nch+1:end) = 0; %% self == 0
F_granger = F

signal_ = zscore(signal,0,2);
F = zeros(nch);
for i = 1:nch
  si = signal_(i,:);
  for j = i+1:nch
    sj = signal_(j,:);
    [F(i,j) F(j,i)] = deal(mutualinfo(si,sj));
  end
end
F_mi = F;
%}

F = F_mi;

[dd,deg_rank] = sort(sum(F));
aa = (averMeas(F));
mm = (moduMeas(F));

chs = cell(8,1);
for i = 1:4, chs{i} = sprintf('CA3(%d)',i); end
for i = 5:8, chs{i} = sprintf('CA1(%d)',i); end

h = plot(deg_rank, aa, 'r.', ...
         deg_rank, mm, 'b.', ...
         'MarkerSize',30, 'LineWidth',3);
set(gca, 'FontSize', 16, ...
         'XLim', [0 8]+.5, ...
         'XTick', 1:8, ...
         'XTickLabel', chs(deg_rank), ...
         'LineWidth',1)
xlabel('Rank of Granger causality')
ylabel('Controllability')
title('Rank of Granger Causality vs Average Controllability')
legend(h, 'Average controllability', 'Modal controllability', 'Location','Best')
% print('-dpng','-r200','results/rank_granger_controllability.png')
