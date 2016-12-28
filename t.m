clear
addpath(genpath('~/src/chronux_2_11/spectral_analysis/'))
addpath(genpath('experiments'))

data = load('data/ARN052_fad_grid_2016_12_20_1_model_data.mat');
biomarker = 'psd';
s_window = [1 50];
window = 0.5;
offset = 0.1;
[xx yy] = gather_training(data, biomarker, s_window, window, offset);

w

% 1. when/where to look to see an effect?
%  - psd over which frequencies?  (s_window)
%  - psd over what time window? (window)
%  - grab window how soon after stimulation? (offset)

% 2. which stimulation parameters have most effect?

return



signal = loadsome('/Volumes/TDT Data/extracted_data/ARN045/4D_single_channel/channel_grid_1_1.mat','data');
signal(:,:,72) = 0;
fs = 24414;
for i = 2:72
  signal(:,:,i) = loadsome(sprintf('/Volumes/TDT Data/extracted_data/ARN045/4D_single_channel/channel_grid_1_%d.mat',i),'data');
end

return

clear
[signal fs stim] = loadcached('data/ARN045_baseline_stim_4V_2000hz','signal','fs','stim');

trial = 3;
channel = 1;

band = [0 500];
sig = signal(:,:,trial);
sig = zscore(sig,0,2);
stim(trial,:)

subplot(2,1,1)
tt = 0:1/fs:(length(sig)-1)/fs;
plot(tt,sig(channel,:)), axis tight, grid on
ylabel('signal')
set(gca, 'FontSize',30);

subplot(2,1,2)
params = struct('Fs', fs, 'tapers', [3 5], 'fpass', band);
[pxx tt ff] = mtspecgramc(sig',[.1 .001],params);
imagesc(tt,ff,10*log10(pxx(:,:,channel))')
axis xy
axis tight
grid on
ylabel('frequency (Hz)')
set(gca, 'FontSize',30);

return

exp_power
return
exp_trials
return

fs = 2000;

[signal stim] = loadcached('data/ARN045_baseline_stim_4V_2000hz','signal','stim');
idx = find(stim.stimulation_channel_order==7);

params = struct('Fs', fs, 'tapers', [6 11], 'fpass', [4 10]);
win = [0.1 0.05]; % window, step (seconds)

lbls = {};
for i = 1:3
  lbls{i} = sprintf('ch%d', i);
end
lbls{4} = sprintf('ch%d', 7);

clf
for i = 1:numel(idx)
  i
  sp(3,1,i)
  % [pxx tt ff] = mtspecgramc(signal([1:3 7],:,idx(i))', win, params);
  % h = plot(tt, 10*log10(pxx));
  h = plot(signal([1:3 7],:,idx(i)));
  legend(h, lbls{:}, 'Location','Best')
  axis tight, axis xy, grid on
  ylabel(sprintf('%d Hz\n', stim.pulse_frequency(idx(i))))
end
return


params = struct('Fs', fs, 'tapers', [6 11], 'fpass', [0 500]);
win = [0.1 0.05]; % window, step (seconds)

pxx = [];
for i = 1:ntrials
  i
  [pxx(:,:,:,i) tt ff] = mtspecgramc(signal(:,:,i)', win, params);
end

pxx_ = mean(pxx,4);
for i = 1:8
  sp(4,2,i)
  imagesc(tt, ff, 10*log10(pxx_(:,:,i))')
  axis tight, axis xy, grid on
  ylabel(sprintf('ch%\n', i))
end


return

clf
exp_predictor_pvalue(data, {'Sub' 'ACg' 'CA1'})
print('-dpng','-r200','results/predictor_pvalue_pos.png')

clf
exp_predictor_pvalue(data, {'TC' 'DLPFC' 'CA3'})
print('-dpng','-r200','results/predictor_pvalue_neg.png')

return

save('results/dloc_pred_linear','dloc')
exp_plot_location_regression(data,'TC',500)
return
% exp_plot_location('EC',250)

return

fs = 6103.52; % sample rate
c_blue = [180 211 248] / 255;

ch = 1;

wn = [100 600]/(fs/2);
[b a] = butter(2, wn);
% sig_ = filtfilt(b,a,sig(:,ch));


params = struct('Fs', fs, 'tapers', [6 11]);
win = [0.05 0.01]; % window, step (seconds)

sig = sig(:,ch);

[pxx tt ff]  = mtspecgramc(sig, win, params);

% imagesc(tt,ff, 10*log10(abs(pxx)'+eps)), axis xy, grid on

mu = ((pxx * ff') ./ sum(pxx,2))';
ff2 = ff.^2;
sd = (sqrt( ((pxx * (ff.^2)') ./ sum(pxx,2)) - (mu').^2 ))';

clf, hold on
h = fill([tt fliplr(tt)], [mu-sd fliplr(mu+sd)], c_blue);
alpha(h, 0.75)
set(h, 'EdgeColor', 'none')

plot(tt,mu,'b')
axis tight, grid on

return

d = data(data.Location=='TC',:);
corrplot([d.Biomarker d.Frequency d.Amplitude d.Duration d.dxBiomarker], ...
         'VarNames',{'Bio' 'Freq' 'Amp' 'Dur' 'dx'})


for ch = 1:8
  sig_pre  = sig(1:fix(30.5*fs),ch);
  sig_post = sig(fix(1200*fs):end,ch);

  % [pxx_pre  tt_pre ff]  = mtspecgramc(sig_pre,  win, params);
  % [pxx_post tt_post ff] = mtspecgramc(sig_post, win, params);

  [pxx_pre  ff_pre ] = mtspectrumc(sig_pre,  params);
  [pxx_post ff_post] = mtspectrumc(sig_post, params);
  pxx_post_ = interp1(ff_post, pxx_post, ff_pre);

  sp(4,2,ch)
  % plot(ff_pre,pxx_pre,'b',  ff_post,pxx_post,'r')
  plot(ff_pre, (pxx_post_ - pxx_pre') ./ abs(pxx_pre'))
  axis tight, grid on
  ylabel(sprintf('ch%d',ch))
  drawnow
end

% c_pre  = mscohere(s6pre,  s11pre,  [], [], 0:150, fs);

% ch = 11;
% sig_pre  = sig(1:fix(30.5*fs),ch);
% sig_post = sig(fix(1200*fs):end,ch);

% c_post = mscohere(s6post, s11post, [], [], 0:150, fs);

return
% sp(2,1,1)
% imagesc(tt_pre,ff,10*log10(abs(pxx_pre)'+eps)); axis xy, grid on
mu_pre = (pxx_pre * ff') ./ sum(pxx_pre,2);
ff2 = ff.^2;
sd_pre = sqrt((pxx_pre * ff2') ./ sum(pxx_pre,2)  - mu_pre.^2);

% sp(2,1,2)
% imagesc(tt_post,ff,10*log10(abs(pxx_post)'+eps)); axis xy, grid on
% xlim([0 max(tt_pre)])
mu_post = (pxx_post * ff') ./ sum(pxx_post,2);
sd_post = sqrt((pxx_post * ff2') ./ sum(pxx_post,2)  - mu_post.^2);

plot(tt_pre, mu_pre, 'b', ...
     tt_post, mu_post, 'r')
legend('pre','post', 'Location','Best')

return

nsamples = numel(lbls);
nchannels = 4; %size(segs{1},2);

wn = [100 600]/(fs/2);
[b a] = butter(2, wn);

sig  = segs{2};
sig = bsxfun(@minus, sig, mean(sig,1));
sig_ = filtfilt(b,a,sig);
assert(~any(isnan(sig_(:))))

mu = mean(sig_,1);
sd = std(sig_,0,1);

tt = 0:1/fs:(size(sig,1)-1)/fs;

for i = 1:nchannels
  sp(2,2,i)
  plot(tt,sig(:,i),'r:', ...
       tt,sig_(:,i),'b', ...
       tt([1 1; end end]), mu(i)+sd(i)*[1 -1; 1 -1]*5, 'y--')
  axis tight, grid on
  ylabel(sprintf('ch%d',i))
end


return


% mtspecgramc parameters
params = struct('Fs', fs, 'fpass', [100 600], 'tapers', [6 11]);
win = [0.2 0.05]; % window, step (seconds)

segs = segs([1 8]);

for i = 1:2

  sig = segs{i};
  sig = bsxfun(@minus, sig, mean(sig,1));
  [Pxx tt ff] = mtspecgramc(sig, win, params);
  isripple = (ff < 250);

  for j = 1:nchannels
    sp(4,2,2*(j-1) + i), hold on
    pxx = Pxx(:,:,j);

    mu = ((pxx * ff') ./ sum(pxx,2))';
    st = (sqrt( ((pxx * (ff.^2)') ./ sum(pxx,2)) - (mu').^2 ))';

    h = fill([tt fliplr(tt)], [mu-st fliplr(mu+st)], c_blue);
    alpha(h, 0.75)
    set(h, 'EdgeColor', 'none')

    plot(tt,mu,'b')
    axis tight, grid on
    ylim([0 600])

    % pxx = 10*log10(abs(pxx)+eps);
    % pxx = bsxfun(@rdivide, pxx, sum(pxx,2));
    % imagesc(tt,ff,pxx');
    % axis tight, grid on

    ylabel(sprintf('#%d, ch%d',i,j))
  end
end

return

pxx_non = cell2mat(psd(:,:,:,lbls==0));
pxx_vol = cell2mat(psd(:,:,:,lbls==1));

for ch = 1:nchannels
  % non-volatile
  sp(nchannels,2,2*(ch-1)+1,0.1)
  mu = mean(pxx_non(:,:,ch,:),4);
  % mu = pxx_non(:,:,ch,2);
  imagesc(tt,ff,10*log10(abs(mu)'+eps))
  colorbar, grid on, axis xy
  ylabel(sprintf('ch%d, Hz',ch))
  xlabel('non-volatile, time (s)')

  % volatile
  sp(nchannels,2,2*(ch-1)+2,0.1)
  mu = mean(pxx_vol(:,:,ch,:),4);
  % mu = pxx_vol(:,:,ch,1);
  imagesc(tt,ff,10*log10(abs(mu)'+eps))
  colorbar, grid on, axis xy
  ylabel(sprintf('ch%d, Hz',ch))
  xlabel('volatile, time (s)')
end

return

feat = {};
parfor i = 1:nsamples
  cx = [];
  for ch1 = 1:nchannels
    for ch2 = ch1+1:nchannels
      cx = [cx; ccx{i,ch1,ch2}];
    end
  end
  feat{i,1} = cx';
end
feat = cell2mat(feat);

for ch1 = 1:nchannels
  for ch2 = ch1+1:nchannels

    vol = cell2mat(ccx(lbls==1,ch1,ch2));
    non = cell2mat(ccx(lbls==0,ch1,ch2));
    mu_vol = mean(vol,1);  st_vol = std(vol,0,1);
    mu_non = mean(non,1);  st_non = std(non,0,1);

    ll = lookup(ch1,ch2);
    sp(3,2,ll,0.07)
    hold on

    [mu st] = deal(mu_non, st_non);
    h = fill([freqs fliplr(freqs)], [mu-st fliplr(mu+st)], c_blue);
    alpha(h, 0.75)
    set(h, 'EdgeColor', 'none')

    [mu st] = deal(mu_vol, st_vol);
    h = fill([freqs fliplr(freqs)], [mu-st fliplr(mu+st)], c_red);
    alpha(h, 0.5)
    set(h, 'EdgeColor', 'none')

    plot(freqs,mu_non,'b', ...
         freqs,mu_vol,'r')

    axis tight, grid on
    ylabel(sprintf('ch%d,ch%d', ch1, ch2))
    xlabel('freqs (Hz)')
  end
end
% print('-dpng','-r300','/tmp/crosscorr')

return


% fs = 6103.52; % sample rate
% c_red = [255 143 146] / 255;
% c_blue = [180 211 248] / 255;

% isvalid = ~isnan(y_hat{1});

% tt = 0:1/fs:(size(segs{1},1)-1)/fs;
% tt = tt(isvalid);
% y_vol = cell2mat(y_hat(lbls==1));  y_vol = y_vol(isvalid,:)';
% y_non = cell2mat(y_hat(lbls==0));  y_non = y_non(isvalid,:)';
% mu_vol = mean(y_vol);   st_vol = std(y_vol);
% mu_non = mean(y_non);   st_non = std(y_non);

% c_red = [255 143 146] / 255;
% c_blue = [180 211 248] / 255;

% hold on

% h = fill([tt fliplr(tt)], [mu_non - st_non fliplr(mu_non+st_non)], c_blue);
% alpha(h, 0.7)
% set(h, 'EdgeColor', c_blue)

% h = fill([tt fliplr(tt)], [mu_vol - st_vol fliplr(mu_vol+st_vol)], c_red);
% alpha(h, 0.5)
% set(h, 'EdgeColor', c_red)

% % plot means last so they show up most
% h = plot(tt, mu_non, 'b', ...
%          tt, mu_vol, 'r', 'LineWidth',2);
% legend(h, 'nonvolatile', 'volatile', 'Location','Best')
% axis tight
% set(gca, 'FontSize',16)
% xlabel('time (s)')
% ylabel('y_{hat}')
% return

% [B fitinfo] = lasso(feat,lbls, 'CV',10, 'Options',statset('UseParallel',true));
b = B(:,fitinfo.Index1SE); % grab recommended coefficients

win = 1; % window (seconds)
dt = 0.25; % step (seconds)
nlags = fix(0.5*fs);

y_hat = cell(size(segs));
parfor j = 1:numel(segs)
  fprintf('worker %d on sample %d\n', get(getCurrentTask(),'ID'), j)
  seg = segs{j}(:,1:4);

  tail = fix(-win*fs):0; % for easier indexing inside loop

  yy = nan(size(seg,1),1);
  for i = fix(win*fs)+1:fix(dt*fs):size(seg,1)
    s = seg(tail+i,:);
    s_win = s(end-nlags:end,:);
    acx = wincorr(s,s_win,nlags);
    x = acx(:);
    for ch1 = 1:4
      for ch2 = ch1+1:4
        cx = wincorr(s(:,ch1),s_win(:,ch2),nlags);
        x = [x; cx];
      end
    end
    yy(i) = b' * x;
  end

  y_hat{j} = yy;
end

tt = 0:1/fs:(size(segs{1},1)-1)/fs;
plot(tt,y_hat{1},'r', tt,y_hat{2},'b')
axis tight, grid on
return

fs = 6103.52; % sample rate

[segs lbls feat] = loadcached('data/ARN038_fpass10_taps611_win1_alpha09','segs','lbls','feat');
[B fitinfo] = lasso(feat,lbls, 'CV',10, 'Alpha',0.9, 'Options',statset('UseParallel',true));
B = B(:,fitinfo.Index1SE); % grab recommended coefficients

if isempty(gcp('nocreate')),  parpool('local'); end

win = 1; % window (seconds)
dt = 0.01; % step (seconds)
params.Fs = fs; % sampling frequency
params.fpass = [0 10]; % min/max frequency
params.tapers = [6 11]; % 3hz bandwidth

y_hat = cell(size(segs));
parfor j = 1:numel(segs)
  fprintf('worker %d on sample %d\n', get(getCurrentTask(),'ID'), j)
  seg = segs{j};

  tail = fix(-win*fs):-1; % for easier indexing inside loop

  yy = nan(size(seg,1),1);
  for i = fix(win*fs)+1:fix(dt*fs):size(seg,1)
    s = seg(tail+i,:);

    s = bsxfun(@minus, s, mean(s,1)); % mean center
                                      % s = locdetrend(s, fs);
    x = mtspectrumc(s, params);

    yy(i) = B' * x(:);
  end

  y_hat{j} = yy;
end

return

arn = '38_3';
if isequal(arn,'38_1')
  fs = 6103.52; % sample rate
  lbls  = loadcached('data/ARN038_pre_stimulation_dataset/ad_labels_038_pre', 'ad_labels_038')';
  lbls(48) = []; % missing stimulation?
  fntemplate = @(id) ...
      sprintf('data/ARN038_pre_stimulation_dataset/ARN038_prestimulation_10s_segment-%d_new',id);
elseif isequal(arn, '35')
  fs = 24414.0625; % sample rate
  lbls = loadcached('data/ARN035_pre_stimulation_dataset/ad_labels_035_pre_new', 'ad_labels_035');
  cycles = lbls(:,1);
  lbls = lbls(:,2);
  fntemplate = @(id) ...
      sprintf('data/ARN035_pre_stimulation_dataset/ARN035_RapidKindling_Block-%d_prestimulation_10s_new',cycles(id));
elseif isequal(arn,'38_2')
  fs = 6103.52; % sample rate
  lbls = loadcached('data/ARN038_pre_stimulation_dataset_2/ad_labeles_038_2', 'ad_labels_038_2');
  fntemplate = @(id) ...
      sprintf('data/ARN038_pre_stimulation_dataset_2/ARN038_prestimulation_10s_segment-%d',id);
elseif isequal(arn,'38_3')
  fs = 6103.52; % sample rate
  lbls = loadcached('data/ARN038_pre_stimulation_dataset_3/ad_labeles_038_3', 'ad_labels_038_3');
  fntemplate = @(id) ...
      sprintf('data/ARN038_pre_stimulation_dataset_3/ARN038_prestimulation_10s_segment-%d',id);
else
  error('invalid arn #%d', arn)
end

nchannels = 8;
nsamples = length(lbls);

% scales = helperCWTTimeFreqVector(1,80,centfrq('morl'),1/fs,4);
% scales = 1:6;

if isempty(gcp('nocreate'))
  parpool('local')
end

feat = {};
segs = {};
parfor i = 1:nsamples
  fprintf('worker %d on sample %d\n', get(getCurrentTask(),'ID'), i)
  fn = feval(fntemplate, i);
  data = double(loadcached(fn, 'data'));
  % feat{i,1} = get_MT_frequency_spectrum_function(data, fs);
  % wave = {};
  % for j = 1:nchannels
  %   wave{j} = sum(abs(cwt(data(:,j), scales, 'morl', 1/fs)),2)';
  % end
  % feat{i,1} = cell2mat(wave);
  segs{i,1} = data;
end
feat = cell2mat(feat);

return

return

psd_ = reshape(psd, nsamples, nbins, nchannels);
for i = 1:nchannels
  sp(4,2,i,0.04)
  ispos = (lbl(1:nsamples) == 1);
  plot(1:nbins, psd_(~ispos,:,i), 'r', ...
       1:nbins, psd_( ispos,:,i), 'w')
  axis tight, box off
end

return


svm_params = {'autoscale', true, ...
              'BoxConstraint', 0.3, ...
              'kernel_function', 'rbf', 'rbf_sigma', 1};
svm = svmtrain(psd_load, lbl, svm_params{:});

idx = -fliplr(1:60381); % to grab 10 sec of data
% nwindows = fix(size(data,2) / fs); % window every 1 second
% ppp = nan(nwindows, nbins*nchannels);
% ccc = nan(nwindows,1);

cycle = 3;
stims = [  915528 ...
          2752687 ...
          4583743 ...
          6414799  ...
          8245855  ...
          10076911 ...
          11907967];

idx_ = stims(cycle) + idx - 2880;
idx_(1)
data_ = data(:,idx_)';
psd_calc = get_MT_frequency_spectrum_function(data_, fs);
% svmclassify(svm, psd)

psd_calc = reshape(psd_calc,      nbins, nchannels);
psd_load = reshape(psd_load,  [], nbins, nchannels);
xx = 1:nbins; %60381;
for i = 1:nchannels
  sp(4,2,i,0.04)
  h = plotyy(1:nbins, [psd_calc(:,i) psd_load(cycle,:,i)'], ...
             1:60381, data_(:,i));
  axis tight, box off
  if i == 1
    legend(h(1), 'calculated', 'stored')
  end
end

return

for i = 1:nwindows
  idx_ = idx + fix((i-0)*fs);
  data_ = data(:,idx_)';
  psd = get_MT_frequency_spectrum_function(data_, fs);

  ppp(i,:) = psd;
  ccc(i) = svmclassify(svm, psd);
  fprintf('T(%.1f%%) data(%d)  %g\n', 100*i/nwindows, idx_(1), ccc(i))

  % for j = 1:2 %nchannels
  %   sp(nchannels,2,2*(j-1)+1,0.04)
  %   plot(data_(j,:))
  %   axis tight, box off
  %   if j == 1
  %     title(T)
  %   end
  %   sp(nchannels,2,2*(j-1)+2,0.04)
  %   plot(psd(:,j))
  %   axis tight, box off
  % end
  % drawnow
end
return






nsamples = numel(lbl);
nch = 8; % number of channels


% svm_params = {'Standardize', true, ...
%               'KernelFunction', 'RBF', ...
%               'KernelScale', 'auto'};
svm_params = {'autoscale', true, ...
              'BoxConstraint', 0.3, ...
              'kernel_function', 'rbf', 'rbf_sigma', 1};
% svm_params = {};

iscorrect = nan(size(lbl));
for i = 1:nsamples
  % prepare new training set by leave-one-out
  [ft_ lbl_] = deal(ft,lbl);
  ft_(i,:) = [];
  lbl_(i)  = [];

  % svm = fitcsvm(ft_, lbl_, svm_params{:});
  svm = svmtrain(ft_, lbl_, svm_params{:});
  iscorrect(i) = (lbl(i) == svmclassify(svm, ft(i,:)));
end
success = mean(iscorrect)

return

% svm_wave = fitcsvm(ft_wave, lbl, svm_params{:});
% cv = crossval(svm_wave,'kfold',5);
% kfoldLoss(cv)

return


clear, clf
% psd  = loadcached('data/ARN038_pre_stimulation_dataset/ad_psd_features_038_pre', 'ad_psd_features_038');
psd  = loadcached('data/ARN038_post_stimulation_dataset/ad_psd_features_038_post', 'ad_psd_features');
lbl  = loadcached('data/ARN038_pre_stimulation_dataset/ad_labels_038_pre', 'ad_labels_038');

c_red = [255 143 146] / 255;
c_blue = [180 211 248] / 255;

nch = 8;
bins = 1:size(psd,2)/nch;
nsamples = size(psd,1);

psd = bsxfun(@minus, psd, mean(psd,2)); % mean shift
psd = bsxfun(@rdivide, psd, sqrt(sum(psd.^2,2))); % normalize

% determine std across all samples to be used in kernel
% K = @(u)  3/4 * (1 - min(1,sqrt(sum(u.^2,2))).^2); % Epanechnikov
K = @(u)  exp(-sum(u.^2,2)/2)/sqrt(2*pi); % Gauss
dd_pos = [];
dd_neg = [];
for i = 1:nsamples
  for j = i+1:nsamples
    if lbl(i) ~= lbl(j)   continue, end

    d = K(psd(i,:) - psd(j,:));
    if lbl(i)==1
      dd_pos(end+1) = d;
    else
      dd_neg(end+1) = d;
    end
  end
end
npos = nnz(lbl==1);
nneg = nnz(lbl==0);
h_pos = (4*std(dd_pos)^5/3/npos)^(1/5)
h_neg = (4*std(dd_neg)^5/3/nneg)^(1/5)

% kern = @(x,y)  exp(-dist(x,y)^2/(2*sigma^2));
kern = @(d,h)  K(d/(2*h));
kern = @(x,y,h)  dot(x,y);

pdf = zeros(2,numel(lbl));
for i = 1:numel(lbl)
  x = psd(i,:);

  pdf_pos = [];
  for j = find(lbl==1)
    if i == j, continue, end % ignore self
    pdf_pos(end+1) = kern(x, psd(j,:), h_pos);
  end
  pdf(1,i) = mean(pdf_pos);

  pdf_neg = [];
  for j = find(lbl==0)
    if i == j, continue, end % ignore self
    pdf_neg(end+1) = kern(x, psd(j,:), h_neg);
  end
  pdf(2,i) = mean(pdf_neg);

end

nn = 1:12;
[lbl(nn); pdf(:,nn); pdf(1,nn)>pdf(2,nn)]

success_rate = 1 - mean(xor(lbl, pdf(1,:)>pdf(2,:)))

return

% something is off with this PRE-stim sample
% psd(44,:) = [];
% lbl(44)   = [];

npos = sum(lbl);
nneg = numel(lbl) - npos;
bins = 1:size(psd,2)/nch;

for ch = 1:nch
  ch_idx = (ch-1)*numel(bins) + bins;
  psd_ch = psd(:,ch_idx);
  psd_pos = psd_ch(lbl==1,:); % after-discharge present
  psd_neg = psd_ch(lbl==0,:); % no after-discharge

  sp(4,2,ch, 0.04)
  hold on

  % positive
  mu_pos = mean(psd_pos);
  st = std(psd_pos);
  % h = fill([bins fliplr(bins)], [mu_pos - st fliplr(mu_pos+st)], c_blue); % 1 std bound
  % alpha(h, 0.5)
  % set(h, 'EdgeColor', 'none')
  plot(bins, psd_pos, 'b:')

  % negative
  mu_neg = mean(psd_neg);
  st = std(psd_neg);
  % h = fill([bins fliplr(bins)], [mu_neg - st fliplr(mu_neg+st)], c_red); % 1 std bound
  % alpha(h, 0.5)
  % set(h, 'EdgeColor', 'none')
  plot(bins, psd_neg, 'r:')

  % plot means last so they show up most
  h = plot(bins, mu_pos, 'b', ...
           bins, mu_neg, 'r', 'LineWidth',2);

  plot(bins, psd_ch(44,:), 'k', 'LineWidth',2)

  axis tight, box off
  xlim([1 16])
  max_mu = max([mu_pos(:); mu_neg(:)]);
  ylim([0 1.5*max_mu])
  ylabel(sprintf('ch%d', ch))

  if ch == 1
    h = legend(h, 'after-discharge present', 'no after-discharge', 'Location', 'NorthEast');
  end
end
return


% addpath(genpath('~/src/chronux_2_11/spectral_analysis'))

data_pre = loadcached('data/ARN038_pre_stimulation_dataset/ARN038_prestimulation_10s_segment-1', 'data');
psd_pre  = loadcached('data/ARN038_pre_stimulation_dataset/ad_psd_features_038_pre', 'ad_psd_features_038');

return

data_post = loadcached('data/ARN038_post_stimulation_dataset/ARN038_poststimulation_10s_segment-1', 'data');
fs = 6103.52; % sample rate

clf

nchannels = size(data_pre,1);
for i = 1:nchannels
  dt = 1/fs;

  time = 0:dt:(size(data_pre,2)-1)*dt;
  sp(nchannels,2,2*(i-1)+1)
  plot(time, data_pre(i,:))
  ylabel(sprintf('ch%d', i))
  axis tight; grid

  time = 0:dt:(size(data_post,2)-1)*dt;
  sp(nchannels,2,2*(i-1)+2)
  plot(time, data_post(i,:))
  ylabel(sprintf('ch%d', i))
  axis tight; grid
end
return

dt = 1/fs;
time = 0:dt:(size(data,2)-1)*dt;

% min/max frequencies
waves.delta = [ 0   4];
waves.theta = [ 4   8];
waves.alpha = [ 8  12];
waves.beta  = [12  30];

% x = data(1,:);
% P_delta = bandpower(x, fs, waves.delta)
% P_theta = bandpower(x, fs, waves.theta)
% P_alpha = bandpower(x, fs, waves.alpha)
% P_beta  = bandpower(x, fs, waves.beta )
% P_total = norm(x,2)^2
% return

nchannels = size(data,1);
for i = 1:nchannels
  eeg = data(i,:);

  sp(nchannels,2,2*(i-1)+1)
  plot(time, eeg)
  ylabel(sprintf('ch%d', i))
  axis tight; grid

  % % construct bandpass filter
  % N = 4; % 4th order filter
  % [c d] = butter(N,(Wn+1)/fs); % construct filter
  % theta = filter(c,d,eeg); % bandpass

  sp(nchannels,2,2*(i-1)+2)
  [pxx f] = pwelch(eeg, [], [], [], fs);
  plot(f, 10*log10(pxx),'y')
  axis tight; grid
end





lbl  = loadcached('data/ARN038_pre_stimulation_dataset/ad_labels_038_pre', 'ad_labels_038') == 1;
nsamples = numel(lbl);
scales = 32;  % scale of wavelet
nch = 8; % number of channels

for i = 1:nsamples
  i
  data = loadcached(sprintf('data/ARN038_pre_stimulation_dataset/ARN038_prestimulation_10s_segment-%d',i), 'data');
  data = double(data); % upcast
  data = bsxfun(@minus, data, mean(data,2)); % mean-shift each channel

  % each channel independently
  ft = [];
  for j = 1:nch
    w = cwt(data(j,:), 1:scales, 'morl', 1/fs);
    ft = [ft; sum(w,2)]; % features are sum of coef weights over the window
  end
  features(i,:) = ft;
end
return


addpath(genpath('~/src/chronux_2_11/spectral_analysis/'))
addpath(genpath('.'))
set(gcf, 'Menubar', 'figure')



lbls = loadcached('data/ARN038_pre_stimulation_dataset_3/ad_labeles_038_3', 'ad_labels_038_3');
fntemplate = @(id) ...
    sprintf('data/ARN038_pre_stimulation_dataset_3/ARN038_prestimulation_10s_segment-%d',id);
for i = 1:numel(lbls)
  fn = fntemplate(i);
  data = loadsome(fn, 'data')';
  save(fn, 'data');
end


feat = feat3;
lbls = lbls3;
npos = nnz(lbls == 1);
nneg = nnz(lbls == 0);

nchannels = 8;
nbins = size(feat,2) / nchannels;
bins = 1:nbins;

c_red = [255 143 146] / 255;
c_blue = [180 211 248] / 255;

for ch = 1:nchannels
  ch_idx = (ch-1)*nbins + bins;
  psd_ch = feat(:,ch_idx);
  psd_pos = psd_ch(lbls==1,:); % after-discharge present
  psd_neg = psd_ch(lbls==0,:); % no after-discharge

  sp(4,2,ch, 0.04)
  hold on

  % positive
  mu_pos = mean(psd_pos);
  st = std(psd_pos);
  h = fill([bins fliplr(bins)], [mu_pos - st fliplr(mu_pos+st)], c_blue); % 1 std bound
  alpha(h, 0.5)
  set(h, 'EdgeColor', 'none')
  % plot(bins, psd_pos, 'b:')

  % negative
  mu_neg = mean(psd_neg);
  st = std(psd_neg);
  h = fill([bins fliplr(bins)], [mu_neg - st fliplr(mu_neg+st)], c_red); % 1 std bound
  alpha(h, 0.5)
  set(h, 'EdgeColor', 'none')
  % plot(bins, psd_neg, 'r:')

  % plot means last so they show up most
  h = plot(bins, mu_pos, 'b', ...
           bins, mu_neg, 'r', 'LineWidth',2);

  axis tight, box off
  xlim([1 16])
  max_mu = max([mu_pos(:); mu_neg(:)]);
  ylim([0 1.5*max_mu])
  ylabel(sprintf('ch%d', ch))

  if ch == 1
    h = legend(h, 'volatile', 'nonvolatile', 'Location', 'NorthEast');
  end
end


return

clf
params = struct('Fs', fs, 'tapers', [6 11]);
win = [0.05 0.01]; % window, step (seconds)

ch = 1;
fs = 2000;
signal = loadcached('data/baseline/control.mat');
s = signal(ch,1:10*fs);

params = struct('Fs', fs, 'fpass', [0 250], 'tapers', [6 11]);
win = [0.2 0.05]; % window, step (seconds)
[pxx tt ff]  = mtspecgramc(s, win, params);

pxx_ = 10*log10(pxx);
imagesc(tt,ff,pxx_');
axis xy
colorbar
grid on
