function exp_plot_location(loc,dur)

  data = loadcached('data/stim_par_table_PS2.1.mat');

  %% Sham - determine mean/stderr as baseline
  d_sham = data(data.Amplitude==0,:);
  xx = table2array(d_sham(:,'dxBiomarker'));
  mu_sham = mean(xx);
  se_sham = stderr(xx);
  nsham = numel(xx);

  %% Stimulation - determine mean/stderr for each frequency/amplitude
  d_stim = data((data.Location == loc) ...
                & (data.Duration == dur | data.Duration == 0) ...
                & (data.Amplitude > 0),:);
  freqs = unique(d_stim.Frequency);
  amps = unique(d_stim.Amplitude);
  nfreqs = numel(freqs);
  namps = numel(amps);
  nsubj = numel(unique(d_stim.Subject));
  [mu se] = deal(nan(nfreqs, namps)); % mean and standard error for each (freq*amp)
  leg = {};
  for i = 1:namps
    d_amp = d_stim(d_stim.Amplitude == amps(i),{'Frequency' 'dxBiomarker'});

    for j = 1:nfreqs
      xx = table2array(d_amp(d_amp.Frequency == freqs(j),'dxBiomarker'));
      mu(j,i) = mean(xx);
      se(j,i) = stderr(xx);
    end
    leg{end+1} = sprintf('%g mA (%d)',amps(i), size(d_amp,1));
  end

  %% plot
  cla, hold on
  % sham
  patch([-10 -10 250 250], ...
        mu_sham + [-se_sham se_sham se_sham -se_sham], ...
        ones(1,3)*0.9, 'EdgeColor', 'none')
  h = errorbar(-5, mu_sham, se_sham);
  % pulse and frequencies
  if nfreqs > 0
    for i = 1:namps    % use loop so errorbar doesn't misunderstand vectors (pulse only)
      offset = (2*i - namps)/2;
      h(end+1) = errorbar(freqs + offset, mu(:,i), se(:,i));
    end
  end
  set(h, 'LineStyle', '--', ...
         'LineWidth', 2, ...
         'Marker', '.', ...
         'MarkerSize', 16);
  set(h(1), 'Color', 'k')

  sham = sprintf('Sham (n=%d)', nsham);
  set(gca, 'XTick', [-5;freqs])
  lbl = get(gca,'XTickLabel');
  lbl(1:2) = {'S' 'P'};
  set(gca,'XTickLabel',lbl, 'FontSize',16)
  plot(xlim, mu_sham([1 1]), 'k--', 'LineWidth', 1)

  ylabel('Change in classifier')
  xlabel('Pulse frequency (Hz)')
  if nfreqs > 0
    set(gca, 'XLim', [-10 210])
  end
  legend(h, sham, leg{:}, 'Location','Best')
  title(sprintf('%s, %d ms, %d subject%s', loc, dur, nsubj, iff(nsubj>1,'s','')))
end
