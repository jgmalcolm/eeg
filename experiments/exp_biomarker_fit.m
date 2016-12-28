clear
data = loadcached('data/stim_par_table_PS2.1.mat');

loc = 'MTL WM';
% subj = 'R1120E'; % !!
% subj = 'R1105E'; % !!!
dur = 500;

%% Stimulation
d_stim = data(  data.Location  ==loc  ...
              & data.Amplitude ==1    ...
              & data.Duration  ==dur,:);

freqs = unique(d_stim.Frequency);
amps = unique(d_stim.Amplitude);
nfreqs = numel(freqs);
namps = numel(amps);

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

%% Sham - determine mean/stderr as baseline
d_sham = data(data.Amplitude == 0,:);
xx = table2array(d_sham(:,'dxBiomarker'));
mu_sham = mean(xx);
se_sham = stderr(xx);
nsham = numel(xx);

%% Plot
clf
hold on
patch([-10 -10 250 250], ...
      mu_sham + 1.96*[-se_sham se_sham se_sham -se_sham], ...
      ones(1,3)*0.9, 'EdgeColor', 'none')
h = errorbar(-5, mu_sham, 1.96*se_sham);
% pulse and frequencies
for i = 1:namps    % use loop so errorbar doesn't misunderstand vectors (pulse only)
  offset = 2*i - namps;
  h(end+1) = errorbar(freqs + offset, mu(:,i), 1.96*se(:,i));
  d_amp = d_stim(d_stim.Amplitude==amps(i),{'Frequency' 'dxBiomarker'});
  hh = plot(d_amp.Frequency, d_amp.dxBiomarker, '.', 'MarkerSize',10);
end
set(h, 'LineStyle', '--', ...
       'LineWidth', 2, ...
       'Marker', '.', ...
       'MarkerSize', 25);
set(h(1), 'Color', 'k')

sham = sprintf('Sham (n=%d)', nsham);
set(gca, 'XTick', [-5;freqs])
lbl = get(gca,'XTickLabel');
lbl{1} = 'S';
if freqs(1) == 0, lbl{2} = 'P'; end
set(gca,'XTickLabel',lbl)
plot(xlim, mu_sham([1 1]), 'k--', 'LineWidth', 1)

ylabel(loc, 'FontWeight','bold')
xlabel('Pulse frequency (Hz)')
if nfreqs > 0
  set(gca, 'XLim', [-10 210])
end


% ignore pulses
d = d_stim(d_stim.Frequency > 0,:);

% only use predictors that have more than one value
preds = {'Frequency' 'Amplitude' 'Duration'};
preds_ = preds(cellfun(@(p) numel(unique(d.(p))) > 1, preds));

dd = d(:,{preds_{:} 'dxBiomarker'});
mdl_lm = fitlm(dd, 'quadratic', ...
            'PredictorVars', preds_, 'ResponseVar','dxBiomarker');

mdl_svm = fitrsvm(dd,'dxBiomarker', ...
              'PredictorNames', preds_, ...
              'KernelFunction','gaussian','KernelScale','auto');

freqs = [10 25 50 100 200];
freqs = [10:10:200];
[y_lm y_lm_ci] = predict(mdl_lm, freqs');
y_svm = predict(mdl_svm, freqs');

hr = plot(freqs,y_svm,'b', 'LineWidth',2);
h(end+1) = errorbar(freqs, y_lm, y_lm_ci(:,1)-y_lm, y_lm_ci(:,2)-y_lm);
set(h, 'LineStyle', '--', ...
       'LineWidth', 2, ...
       'Marker', '.', ...
       'MarkerSize', 25);

legend([h(:);hr(:)], sham, leg{:}, 'quadratic reg', 'SVM reg', 'Location','Best')

set(gca, 'FontSize',14)
return
y = d.dxBiomarker;
% y_hat = predict(mdl, d);
y_hat = resubPredict(mdl);
xx = 1:numel(y);

[~,idx] = sort(y);
plot(xx,y(idx),'b+', xx,y_hat(idx),'r.', xx(mdl.IsSupportVector),y(mdl.IsSupportVector),'mo')

return
%{


for loc = 1:nlocations

  %% Stimulation - determine mean/stderr for each frequency/amplitude
  d_stim = data((data.Location == locations(loc)) ...
                & (data.Duration == dur | data.Duration == 0) ...
                & (data.Amplitude > 0),:);
  freqs = unique(d_stim.Frequency);
  amps = unique(d_stim.Amplitude);
  nfreqs = numel(freqs);
  namps = numel(amps);
  [mu se] = deal(nan(nfreqs, namps)); % mean and standard error for each (freq*amp)
  leg = {};
  for i = 1:namps
    d_amp = d_stim(d_stim.Amplitude == amps(i),{'Frequency' 'dxBiomarker'});

    for j = 1:nfreqs
      xx = table2array(d_amp(d_amp.Frequency == freqs(j),'dxBiomarker'));
      mu(j,i) = mean(xx);
      se(j,i) = stderr(xx);

      % outside sham?
      low  = @(mu,se) mu-1.96*se;
      high = @(mu,se) mu+1.96*se;
      if high(mu_sham,se_sham) < low(mu(j,i),se(j,i)) ...
        || high(mu(j,i),se(j,i)) < low(mu_sham,se_sham)
        fprintf('%4d ms, %3d Hz, %.2f mA, %s\n', ...
                dur, freqs(j), amps(i), char(locations(loc)));
        interesting = [interesting; table(dur, freqs(j), amps(i), locations(loc))];
      end
    end
    leg{end+1} = sprintf('%g mA (%d)',amps(i), size(d_amp,1));
  end

  %% plot
  sp(6,4,loc,0.1)
  hold on
  % sham
  patch([-10 -10 250 250], ...
        mu_sham + 1.96*[-se_sham se_sham se_sham -se_sham], ...
        ones(1,3)*0.9, 'EdgeColor', 'none')
  h = errorbar(-5, mu_sham, 1.96*se_sham);
  % pulse and frequencies
  if nfreqs > 0
    for i = 1:namps    % use loop so errorbar doesn't misunderstand vectors (pulse only)
      offset = 2*i - namps;
      h(end+1) = errorbar(freqs + offset, mu(:,i), 1.96*se(:,i));
    end
  end
  set(h, 'LineStyle', '--', ...
         'LineWidth', 1, ...
         'Marker', '.', ...
         'MarkerSize', 16);
  set(h(1), 'Color', 'k')

  sham = sprintf('Sham (n=%d)', nsham);
  set(gca, 'XTick', [-5;freqs])
  lbl = get(gca,'XTickLabel');
  lbl(1:2) = {'S' 'P'};
  set(gca,'XTickLabel',lbl)
  plot(xlim, mu_sham([1 1]), 'k--', 'LineWidth', 1)

  ylabel(char(locations(loc)), 'FontWeight','bold')
  xlabel('Pulse frequency (Hz)')
  if nfreqs > 0
    set(gca, 'XLim', [-10 210])
  end
  legend(h, sham, leg{:}, 'Location','Best', 'Orientation','horizontal')
  drawnow
end

%}
