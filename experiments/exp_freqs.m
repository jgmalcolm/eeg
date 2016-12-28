clear
[signal fs] = loadcached('data/ARN045_baseline_16.mat','signal','fs');
nch = size(signal,1);

ntrials = 50; % number of trials
maxorder = 25; % maximum model order to consider
acmaxlags = 0; % estimate from data

% recording channels
% lookup = [2 4 6 8 9 11 13 15];
lookup = [9 8 10 7 11 6 12 5 13 4 14 3 15 2 16 1];

wins = 2

nwins = numel(wins);
for i = 1:nwins
  win = wins(i);
  fprintf('#%d %.2f\n', pindex, win)

  trial_len = fix(win * fs); % sec

  strials = [];
  for j = 1:ntrials
    idx = (j-1)*trial_len + (1:trial_len);
    for k = 1:nch
      strials(k,:,j) = detrend(signal(k,idx));
    end
  end

  [GCf,~,~,fres] = est_grangerfreq(strials, maxorder, acmaxlags);

  ff = sfreqs(fres, fs);

  stim16 = [16 12 8 4 3 7 11 15];
  read16 = [14 10 6 2 1 5 9 13];
  stim8  = [1 3 5 7 10 12 14 16];
  read8  = [2 4 6 8 9 11 13 15];
  nsch = numel(stim8);
  clf
  ymax = [];
  for r = 1:nsch
    for c = 1:nsch
      sp(nsch,nsch, (r-1)*nsch+c)
      h = plot(ff, squeeze(GCf(read16(c),stim16(r),:)), ...
               ff, squeeze(GCf(stim16(r),read16(c),:)), 'LineWidth',3);
      set(h(1),'Color',.8*[1 1 1], 'LineStyle',':')
      axis tight
      xlim([50 350])
      grid on, box off
      if r+1==c
        set(gca, 'YTickLabel','')
      else
        set(gca, 'XTickLabel','','YTickLabel','')
      end
      set(gca, 'FontSize',15)
      ymax(end+1) = max(ylim);
    end
  end

  for r = 1:nsch
    for c = 1:nsch
      sp(nsch,nsch, (r-1)*nsch+c)
      set(gca, 'YLim', [0 3*mean(ymax)])
      text(60, .9*max(ylim), sprintf('%d to %d',stim8(c),read8(r)), 'FontSize',20)
    end
  end

  savefig(gcf, 'results/freqs.fig','compact')
  print('-dpng','-r200',sprintf('/tmp/freqs_%.1f.png', win))

end
