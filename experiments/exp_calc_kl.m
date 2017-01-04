function kl = exp_calc_kl(data, windows, offsets)

noffsets = numel(offsets);
nwindows = numel(windows);
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
    fprintf('%s %.0f%%\n', mfilename, ((i-1)*noffsets + j)/(noffsets*nwindows)*100)
  end
end
vars = {'Sample' 'Window' 'Offset' 'Divergence'};
kl.Properties.VariableNames = vars;
