clear
mse = loadcached('results/dloc_global');
n = size(mse,1);

loc = mse.Location;
X = [mse.mseGlobal mse.mseLocal];
% reorder
loc = loc(end:-1:1);
X = X(end:-1:1,:);

h = barh(X);
set(gca, 'YTick', 1:n, ...
         'YTickLabel', char(loc), ...
         'YLim', [0 n]+.5, ...
         'FontSize', 16)
set(h, 'BarWidth', 1, 'EdgeAlpha',0)
set(h(1), 'FaceColor', hex2rgb('418CF0')/256)
set(h(2), 'FaceColor', hex2rgb('FCB441')/256)
legend(h, 'Global', 'Local', 'Location','Best')
xlabel('MSE')
title('Local vs Global')

print('-dpng', '-r200', 'results/local_vs_global_bar.png')
