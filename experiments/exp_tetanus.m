clear
addpath('~/src/Network-Controllability-Diagnostics/')
addpath('~/src/mutualinfo/')
fs = 6103.52;
s_ctl = loadcached('data/baseline/control.mat');
[s_tet1 s_tet2 s_tet3] = loadcached('data/baseline/tetanus.mat','x1','x2','x3');

% s_ctl = zscore(s_ctl,0,2);
% s_tet = zscore(s_tet,0,2);

warning('using reduced size dataset')
tsec = 4;
s_ctl = s_ctl(:,1:fix(tsec*fs));
s_tet = s_tet1(:,1:fix(tsec*fs));

maxorder = 25;
acmaxlags = 1000;
X_ctl = est_granger(s_ctl, maxorder, acmaxlags);
X_tet = est_granger(s_tet, maxorder, acmaxlags);

% X_ctl = est_mutualinfo(s_ctl);
% X_tet = est_mutualinfo(s_tet3);

inf_ctl = sum(X_ctl);
inf_tet = sum(X_tet);

avg_ctl = averMeas(X_ctl);
avg_tet = averMeas(X_tet);


h = plot(inf_ctl(1:8),  avg_ctl(1:8), 'r.', ...
         inf_ctl(9:16), avg_ctl(9:16), 'b.', ...
         inf_tet(1:8),  avg_tet(1:8), 'k.', ...
         inf_tet(9:16), avg_tet(9:16), 'm.', ...
         'MarkerSize',30, 'LineWidth',3);
set(gca, 'FontSize', 16, ...
         'LineWidth',1)
xlabel('Mutual information')
ylabel('Average controllability')
title('Mutual Information vs Average Controllability')
legend(h, 'Control CA3', 'Control CA1', 'Tetanus CA3', 'Tetanus CA1', 'Location','Best')

print('-dpng','-r200','results/tetanus_gc_tet1.png')
