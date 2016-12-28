addpath('~/src/Network-Controllability-Diagnostics')
clear
fs = 2000;
signal = loadcached('data/baseline/tetanus.mat','x1');
signal = signal(1:2:end,:); % drop even channels for computational speed
nch = size(signal,1);

fs_ = 1000;
signal_ = [];
for i = 1:nch
  signal_(i,:) = resample(signal(i,:), fs_, fs);
end
signal = signal_;
fs = fs_;

maxorder = 25; % maximum model order to consider
acmaxlags = 0; % estimate from data

winlim = [55 60]; % window to examine (sec)
winlim = [57 58]
win = 0.5; % sec
dt = 0.1; % time step
tt = winlim(1):dt:winlim(2);
nwin = numel(tt)


[mo ac] = deal(nan(nwin,1));
[gc mi] = deal(cell(1,1,nwin));
parfor i = 1:nwin
  ts = winlim(1) + (i-1)*dt;
  idx = fix((winlim(1) + (i-1)*dt)*fs) + (-fix(win*fs):0);
  s = signal(:,idx);
  tic
  [g mo(i) ac(i)] = est_granger(s, maxorder, acmaxlags);
  m = est_mutualinfo(s);
  T = toc;
  gc{i} = g;
  mi{i} = m;
  fprintf([sprintf('#%d:%.2f  order %d  lags %4d  %.0f sec\n  gc', ...
                   get(getCurrentTask(),'ID'), ts, mo(i), ac(i), T) ...
           sprintf('  %.5f', sum(g,1)) ...
           '\n  mi' ...
           sprintf('  %.5f', sum(m,1)) ...
           '\n'])
end

gc_ = cell2mat(gc);
mi_ = cell2mat(mi);

clf

sp(3,1,1)
tt = 0:1/fs:(length(signal)-1)/fs;
idx = fix(winlim(1)*fs):fix(winlim(2)*fs);
h = plot(tt(idx), signal(:,idx), 'LineWidth', 3);
ylabel('Signal')
axis tight, grid on
set(gca, 'FontSize', 22)
lbl = {};
for i = 1:nch
  lbl{i} = sprintf('channel %d (%s)', ...
                   1+2*(i-1), ... % selected only odd channels above
                   iff(i>4,'CA1','CA3'));
end
legend(h,lbl{:}, 'Location','West')
set(h(5:8),'LineStyle','--')

tics = get(gca, 'XTick');
lbls = get(gca, 'XTickLabel');

sp(3,1,2)
tt = winlim(1):dt:winlim(2);
h = plot(tt, squeeze(sum(gc_,1)), 'LineWidth', 3);
xlabel('Time (seconds)')
ylabel('Granger Causality')
axis tight, grid on
set(gca, 'FontSize', 22, 'XTick',tics, 'XTickLabelMode','auto')
set(h(5:8),'LineStyle','--')


% pull together channels
for i = 1:nch, mi_(i,i,:) = 0; end % zero self
sp(3,1,3)
h = plot(tt, squeeze(sum(mi_,1)), 'LineWidth', 3);
axis tight, grid on
xlabel('Time (seconds)')
ylabel('Mutual Information')
set(gca, 'FontSize', 22, 'XTick',tics, 'XTickLabelMode','auto')
set(h(5:8),'LineStyle','--')


savefig(gcf, 'results/gcmi_seizure.fig','compact')
print('-dpng','-r200','/tmp/gcmi_seizure.png')
