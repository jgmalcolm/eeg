clear
data = loadcached('data/stim_par_table_PS2.1.mat');

% Sham - determine mean/stderr as baseline
d_sham = data(data.Amplitude == 0,:);
mu_sham = mean(d_sham.dxBiomarker);
se_sham = stderr(d_sham.dxBiomarker);


% grab stimulations
d_stim = data(data.Amplitude >0,:);

locations = unique(d_stim.Location);   nlocs = numel(locations);
durs = unique(d_stim.Duration);        ndurs = numel(durs);
freqs = unique(d_stim.Frequency);      nfreqs = numel(freqs);
amps = unique(d_stim.Amplitude);       namps = numel(amps);

[ff aa dd] = deal(nan(1,nlocs));
for loc = 1:nlocs
  d_loc = d_stim(d_stim.Location == locations(loc),:);

  is_above = @(vv) (vv - mu_sham)/se_sham >  1.96;
  is_below = @(vv) (vv - mu_sham)/se_sham < -1.96;

  %% frequency
  [num den] = deal(0);
  for i = 1:nfreqs
    vv = d_loc.dxBiomarker(d_loc.Frequency == freqs(i));
    nabove = nnz(is_above(vv));
    nbelow = nnz(is_below(vv));
    num = num + max(nabove, nbelow);
    den = den + numel(vv);
  end
  ff(loc) = num/den;

  %% amplitude
  [num den] = deal(0);
  for i = 1:namps
    vv = d_loc.dxBiomarker(d_loc.Amplitude == amps(i));
    nabove = nnz(is_above(vv));
    nbelow = nnz(is_below(vv));
    num = num + max(nabove, nbelow);
    den = den + numel(vv);
  end
  aa(loc) = num/den;

  %% duration
  [num den] = deal(0);
  for i = 1:ndurs
    vv = d_loc.dxBiomarker(d_loc.Duration == durs(i));
    nabove = nnz(is_above(vv));
    nbelow = nnz(is_below(vv));
    num = num + max(nabove, nbelow);
    den = den + numel(vv);
  end
  dd(loc) = num/den;
end

d = {'Frequency' ff; ...
     'Duration'  dd; ...
     'Amplitude' aa};
for i = 1:3
  lbl = d{i,1};
  xx = d{i,2};

  subplot(1,3,i)
  [~,idx] = sort(xx);
  plot(xx(idx),1:nlocs,'r', 'Marker','.','MarkerSize',16);
  ylim([0 nlocs]+0.5)
  hold on
  plot([0 0], ylim, 'k--')
  hold off
  xlabel(sprintf('corr(%s,dxBio)',lbl))
  title(lbl)
  set(gca, 'YTick', 1:nlocs, ...
           'YTickLabel', char(locations(idx)), ...
           'FontSize', 12)
  grid on
end
