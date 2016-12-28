addpath('~/src/Network-Controllability-Diagnostics')
clear
fs = 2000;
signal = loadcached('data/baseline/control.mat');
signal = signal(1:2:end,:); % drop even channels
nch = size(signal,1);

maxorder = 25; % maximum model order to consider
acmaxlags = 2000; % maximum autocorrelation lags

win = [0.1:.1:3]; % window size (seconds)
nwin = numel(win);

gc = cell(1,1,nwin);
mo = zeros(nwin,1);
parfor i = 1:nwin
  w = win(i);
  s = signal(:,1:fix(w*fs));
  [g mo(i)] = est_granger(s, maxorder, acmaxlags);
  gc{i} = g;

  [w mo(i) sum(g,1)]
end

gc = cell2mat(gc);

clf
h = plot(win, squeeze(sum(gc,1)));
xlabel('Window (seconds)')
ylabel('Total Granger Causality')
for i = 1:nch
  lbl{i} = sprintf('channel %d (%s)', ...
                   1+2*(i-1), ... % selected only odd channels above
                   iff(i>4,'CA1','CA3'));
end
legend(h,lbl{:}, 'Location','Best')
axis tight
set(gca, 'FontSize', 22)
set(h, 'LineWidth', 3)
set(h(5:8),'LineStyle','--')
