function bg = exp_calc_kl_bg(data, windows)

nwindows = numel(windows);
nsamples = length(data.time_data);
offset = 0; % post-stim offset does not matter
bg = table();
for i = 1:nwindows
  window = windows(i);
  [Spre ~] = gather_prepost(data, window, offset);
  [sum_x sum_xx n] = deal(0);
  for sample = 1:nsamples
    for sample2 = 1:nsamples
      if sample == sample2, continue, end
      d = est_divergence(Spre(:,:,sample), Spre(:,:,sample2));
      sum_x = sum_x + d;
      sum_xx = sum_xx + d*d;
      n = n + 1;
    end
  end
  Ex = sum_x / n;
  Exx = sum_xx / n;
  mu = Ex;
  st = sqrt(Exx - Ex^2);
  bg = [bg; table(window, mu, st)]
  fprintf('progress %.0f%%\n', i/nwindows*100)
end
vars = {'Window' 'Mean' 'Std'};
bg.Properties.VariableNames = vars;
