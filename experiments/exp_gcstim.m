clear
[signal stim fs] = loadcached('data/ARN045_baseline_stim_4V_2000hz', 'signal', 'stim', 'fs');
[nch,~,ntrials] = size(signal);

%% extract pre/stim/post segments
sig_pre  = signal(:,1:4000,:);
sig_stim = signal(:,4001:8000,:);
% sig_post = signal(:,8001:end,:);

maxorder = 25; % maximum model order to consider
acmaxlags = 0; % estimate from data
fmax = []; % max frequency for GC

ch = 1;
freqs = unique(stim.Frequency);

for i = 1:numel(freqs)

  freq = freqs(i);
  idx = (stim.Channel==ch)&(stim.Frequency==freq);
  GCf = est_grangerfreq(sig_pre(:,:,idx), maxorder, acmaxlags, fmax);

  % recording channels
  lookup = [2 4 6 8 9 11 13 15];

  clf
  for r = 1:nch
    for c = 1:nch
      if r == c, continue, end
      sp(nch,nch, (r-1)*nch+c)
      plot(squeeze(GCf(r,c,:)), 'LineWidth',2)
      axis tight
      xlim([0 400])
      ylim([0 .5])
      grid on, box off
      if r+1==c
        xlabel('freq (Hz)')
        ylabel('GC')
      else
        set(gca, 'XTickLabel','')
      end
      set(gca, 'YTickLabel','')
      text(50, 0.85*max(ylim), sprintf('%d to %d',lookup(c),lookup(r)), 'FontSize',12)
    end
  end

  print('-dpng','-r200',sprintf('/tmp/gcf_%d.png', freq))
end
