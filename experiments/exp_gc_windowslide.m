addpath('~/src/Network-Controllability-Diagnostics')
clear
fs = 2000;
signal = loadcached('data/baseline/control.mat');
signal = signal(1:2:end,:); % drop even channels
nch = size(signal,1);

maxorder = 25; % maximum model order to consider
acmaxlags = 2000; % maximum autocorrelation lags

win = 5;
widx = 1:fix(win*fs);
dt = 1; % time step
nwin = fix((size(signal,2)-fs*win) / fs / dt);
nwin = 180

mo = nan(nwin,1);
gc = cell(1,1,nwin);
parfor i = 1:nwin
  idx = (i-1)*dt*fs + widx;
  s = signal(:,idx);
  [g mo(i)] = est_granger(s, maxorder, acmaxlags);
  gc{i} = g;
  [mo(i) sum(g,1)]
end

gc_ = cell2mat(gc);

clf
tt = ((1:nwin)-1)*dt;
h = plot(tt, squeeze(sum(gc_,1)));
xlabel('Time (seconds)')
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

savefig(gcf, sprintf('results/gc_windowslide_%ds.fig',win),'compact')
print('-dpng','-r200',sprintf('/tmp/gc_windowslide_%ds.png',win))
