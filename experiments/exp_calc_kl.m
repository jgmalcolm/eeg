% data = loadcached('data/ARN053_fad_grid_model_data');
data = loadcached('data/ARN052_fad_grid_2016_12_20_1_model_data');

offsets = 0:.1:2;   noffsets = numel(offsets);
windows = 0.1:0.1:1.5;  nwindows = numel(windows);
kl = table();
for i = 1:nwindows
  window = windows(i);
  for j = 1:noffsets
    offset = offsets(j);
    [Spre Spost] = gather_prepost(data, window, offset);
    for sample = 1:size(Spre,3)
      div = est_divergence(Spre(:,:,sample), Spost(:,:,sample));
      kl = [kl; table(sample, window, offset, div)];
    end
    fprintf('progress %.0f%%\n', ((i-1)*noffsets + j)/(noffsets*nwindows)*100)
  end
end
vars = {'Sample' 'Window' 'Offset' 'Divergence'};
kl.Properties.VariableNames = vars;
