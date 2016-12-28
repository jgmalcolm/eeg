clear
[lin qua step] = loadcached('results/dloc_predict','dloc_linear','dloc_quadratic','dloc_stepwise');
n = size(lin,1);

calc_mse = @(p1,p2) cellfun(@(y,y_hat) mean((y-y_hat).^2), p1,p2);

locs = lin.Location;
X = [calc_mse(lin.Y, lin.Y_stim) ...
     calc_mse(lin.Y, lin.Y_bio) ...
     calc_mse(lin.Y, lin.Y_biostim)];
% reorder
[X order] = sortrows(X,-3);
locs = locs(order);

h = barh(X);
set(gca, 'YTick', 1:n, ...
         'YTickLabel', char(locs), ...
         'YLim', [0 n]+.5, ...
         'FontSize', 16)
set(h, 'BarWidth', 1, 'EdgeAlpha',0)
set(h(1), 'FaceColor', hex2rgb('418CF0')/256)
set(h(2), 'FaceColor', hex2rgb('FCB441')/256)
set(h(3), 'FaceColor', hex2rgb('2E8B57')/256)
legend(h, 'Stimulation', 'Biomarker', 'Stim+Bio', 'Location','Best')
xlabel('MSE')
title('MSE using Stimulation, Biomarker, or Bio+Stim')

print('-dpng', '-r200', 'results/predict_biostim_bar.png')
