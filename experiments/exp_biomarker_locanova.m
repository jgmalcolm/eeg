clear
data = loadcached('data/stim_par_table_PS2.1.mat');

locations = unique(data.Location);
nlocs = numel(locations);

% Sham - determine mean/stderr as baseline
d_sham = data(data.Amplitude == 0,:);
median_sham = median(d_sham.dxBiomarker);

vars = {'Location' 'Frequency' 'Amplitude' 'Duration' 'Subjects' 'Samples' 'Pvalue'};
[pvals_pos pvals_neg] = deal(table());

for loc = 1:nlocs
  d_stim = data(data.Location==locations(loc) & data.Amplitude>0,:);

  freqs = unique(d_stim.Frequency)'; nfreqs = numel(freqs);
  amps = unique(d_stim.Amplitude)'; namps = numel(amps);
  durs = unique(d_stim.Duration)'; ndurs = numel(durs);

  [pp_pos pp_neg] = deal(cell(0, numel(vars)));
  ii = 1;
  for freq = freqs
    for amp = amps
      for dur = durs
        d = d_stim(d_stim.Frequency==freq ...
                   & d_stim.Amplitude==amp ...
                   & d_stim.Duration==dur, :);
        nsubj = numel(unique(d.Subject));
        nsamples = numel(d.dxBiomarker);

        %% positive
        if nsamples
          p = signrank(d.dxBiomarker, median_sham, 'tail','right');
        else
          p = nan;
        end
        pp_pos(ii,:) = {locations(loc) freq amp dur nsubj nsamples p};

        %% negative
        if nsamples
          p = signrank(d.dxBiomarker, median_sham, 'tail','left');
        else
          p = nan;
        end
        pp_neg(ii,:) = {locations(loc) freq amp dur nsubj nsamples p};

        ii = ii + 1;
      end
    end
  end

  % remove NaN and append rows
  pp_pos(cellfun(@(pval) isnan(pval), pp_pos(:,end)),:) = [];
  pp_neg(cellfun(@(pval) isnan(pval), pp_neg(:,end)),:) = [];
  pvals_pos = [pvals_pos; cell2table(pp_pos)];
  pvals_neg = [pvals_neg; cell2table(pp_neg)];
end

pvals_pos.Properties.VariableNames = vars;
pvals_neg.Properties.VariableNames = vars;

[pp_pos pp_neg] = deal(nan(1,nlocs));
for loc = 1:nlocs
  p = pvals_pos(pvals_pos.Location==locations(loc),:).Pvalue;
  p_ = bonf_holm(p); % correct for multiple comparisons
  p_ = min(1,p_); % BUG overcorrected?
  pp_pos(loc) = prod(p_);

  p = pvals_neg(pvals_neg.Location==locations(loc),:).Pvalue;
  p_ = bonf_holm(p); % correct for multiple comparisons
  p_ = min(1,p_); % BUG overcorrected?
  pp_neg(loc) = prod(p_);
end

%% Plot

subplot(1,2,1)
[~,idx] = sort(pp_neg,'descend');
plot(pp_neg(idx),1:nlocs,'r', 'Marker','.','MarkerSize',20);
ylim([0 nlocs]+0.5)
hold on
plot([0.05 0.05], ylim, 'k--')
hold off
xlabel('min p-value')
title('locations ranked by smallest p-value for negative effect')
set(gca, 'XLim', [0 .1], ...
         'XTick', [0 .01 .05 .1], ...
         'YTick', 1:nlocs, ...
         'YTickLabel', char(locations(idx)), ...
         'FontSize', 16)
grid on

subplot(1,2,2)
[~,idx] = sort(pp_pos,'descend');
plot(pp_pos(idx),1:nlocs,'r', 'Marker','.','MarkerSize',20);
ylim([0 nlocs]+0.5)
hold on
plot([0.05 0.05], ylim, 'k--')
hold off
xlabel('min p-value')
title('locations ranked by smallest p-value for positive effect')
set(gca, 'XLim', [0 .1], ...
         'XTick', [0 .01 .05 .1], ...
         'YTick', 1:nlocs, ...
         'YTickLabel', char(locations(idx)), ...
         'FontSize', 16)
grid on
