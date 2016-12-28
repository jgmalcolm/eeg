clear

[signal stim fs] = loadcached('data/ARN045_baseline_stim_4V_2000hz', 'signal', 'stim', 'fs');
[nch,~,ntrials] = size(signal);

sig_pre  = signal(:,1:fix(1.5*fs),:);
sig_stim = signal(:,fix(2.5*fs):fix(3.5*fs),:);
sig_post = signal(:,fix(4.5*fs):end,:);

freqs = [100 200 300];
for i = 1:numel(freqs)
  freq = freqs(i);

  idx = find(stim.Frequency==freq);
  nidx = numel(idx);
  params = struct('Fs', fs, 'tapers', [3 5], 'fpass', freq+[-5 5]);

  sig = sig_pre(:,:,idx);
  s = reshape(permute(sig,[2 1 3]), [], nch*nidx);
  % pxx = flat(10*log10(mtspectrumc(s, params)));
  pxx = 10*log10(bandpower(s, fs, freq+[-5 5]));
  pxx_pre = pxx;

  sig = sig_stim(:,:,idx);
  s = reshape(permute(sig,[2 1 3]), [], nch*nidx);
  % pxx = flat(10*log10(mtspectrumc(s, params)));
  pxx = 10*log10(bandpower(s, fs, freq+[-5 5]));
  pxx_stim = pxx;

  sig = sig_post(:,:,idx);
  s = reshape(permute(sig,[2 1 3]), [], nch*nidx);
  % pxx = flat(10*log10(mtspectrumc(s, params)));
  pxx = 10*log10(bandpower(s, fs, freq+[-5 5]));
  pxx_post = pxx;

  mu(1,i) = mean(pxx_pre);
  sd(1,i) = std(pxx_pre);

  mu(2,i) = mean(pxx_stim);
  sd(2,i) = std(pxx_stim);

  mu(3,i) = mean(pxx_post);
  sd(3,i) = std(pxx_post);
end

mu,sd

clf
h = errorbar([0 1 2]'*[1 1 1]+[1 1 1]'*[-1 0 1]/25, mu, sd, 'LineWidth',4, 'Marker','.','MarkerSize',40);
set(gca, 'XTick', [0 1 2], 'XTickLabel', {'pre-stim' 'stim' 'post-stim'}, 'FontSize',40)
ylabel('band power (dB)')

legend(h, '100 Hz', '200 Hz', '300 Hz', 'Location','Best')
