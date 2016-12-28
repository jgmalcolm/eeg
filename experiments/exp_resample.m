addpath('~/src/Network-Controllability-Diagnostics')
clear

fs = 2000;
signal = loadcached('data/baseline/control.mat');
signal = signal(1:2:end,:); % drop even channels
nch = size(signal,1);

fs_ = 250;
signal_ = resample(signal', fs_, fs)'; % downsample

maxorder = 25; % maximum model order to consider
acmaxlags = 2000; % maximum autocorrelation lags

win = 3;

s = signal(:,1:win*fs);
tic
[g mo] = est_granger(s, maxorder, acmaxlags);
toc

s_ = signal_(:,1:win*fs_);
tic
[g_ mo_] = est_granger(s_, maxorder, acmaxlags);
toc

mo,mo_
