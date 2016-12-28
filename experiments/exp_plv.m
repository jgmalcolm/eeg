  eegData = rand(28, 1000, 231);
  srate = 500; %Hz
  filtSpec.order = 50;
  filtSpec.range = [35 45]; %Hz
  dataSelectArr = rand(231, 1) >= 0.5; % attend trials
  dataSelectArr(:, 2) = ~dataSelectArr(:, 1); % ignore trials
  [plv] = pn_eegPLV(eegData, srate, filtSpec, dataSelectArr);
  figure; plot((0:size(eegData, 2)-1)/srate, squeeze(plv(:, 17, 20, :)));
  xlabel('Time (s)'); ylabel('Plase Locking Value');
return


fs = 6103.52; % sample rate
sig = loadcached('data/ARN038_Grid_Search_MainDataTank_Block-50','data');
sig = sig(1,1:fix(fs*1))';
sig = bsxfun(@minus, sig, mean(sig));

tt = 0:1/fs:(size(sig)-1)/fs;

% band-pass
frange = [35 45]; % Hz (gamma band)

order = 4;
wn = frange/(fs/2);
[b a] = butter(order, wn);
sig_butter = filtfilt(b,a,sig);
assert(~any(isnan(sig_butter)))

order = 1000;
f = fir1(order, 2/fs*frange);
sig_fir1 = filter(f, 1, sig);
assert(~any(isnan(sig_fir1)))

subplot(2,1,1)
plot(tt,sig), axis tight, grid on

subplot(2,1,2)
plot(tt,sig_butter,'r',  tt,sig_fir1,'b'), axis tight, grid on

return
subplot(4,1,3)
h_ = hilbert(sig_);
plot(tt,real(h_),'r',  tt,imag(h_),'b'), axis tight, grid on

axis tight
