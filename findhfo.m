function [slow fast] = findhfo(sig, fs)
% Find indices of High Frequency Oscillations
%
% [SLOW FAST] = FINDHFO(SIG,FS) returns indices for slow and fast
% ripples in signal SIG sampled at frequency FS.  SIG must be column
% vector.
%
% - band-pass (100-600 Hz) filtered using a FIR filter.
% - calculate standard deviation across entire signal.
% - sliding window of 200ms moved along the signal in steps of
%   200ms.
% - for each window with std(win) > 2*std(sig), select spikes with
%   amplitude 5x*std(sig)
% - if more than 6 spikes detected, calculate median interspike
%   interval
% - if interval < 4ms then 'fast ripple', if 4-10ms then 'slow ripple'
%
% See Jiruska et al. "Epileptic high-frequency network activity in a
% model of non-lesional temporal lobe epilepsy" in Brain
% 133:1380-1390, 2010.

  assert(iscolumn(sig))

  % parameters from Jiruska2010
  param.frange = [100 600]; % frequency range to examine (Hz)
  param.win = 0.200; % sliding window (sec)
  param.interesting = 1; % std multiple for window to be interesting
  param.amplitude = 1; % amplitude for spikes to be interesting
  param.nspikes = 6; % number of spikes to be detectable
  param.slow = [0.004 0.010]; % interval to be considered slow ripple (sec)
  param.fast = [0.000 0.004]; % interval to be considered fast ripple (sec)

  % band-pass
  wn = param.frange/(fs/2);
  [b a] = butter(2, wn);
  sig = filtfilt(b,a,sig);
  assert(~any(isnan(sig)))

  % examine each window
  sd_sig = std(sig);
  step = fix(param.win*fs);
  [slow fast] = deal([]);
  for i = 1:step:numel(sig)
    swin = sig(i:min(i+step-1,end));
    sd_win = std(swin);

    % interesting window?
    if sd_win <= param.interesting * sd_sig, continue, end

    % find big spikes
    is_big = abs(swin) > param.amplitude * sd_sig;
    dx = diff(swin,1);
    is_signchange = [false; sign(dx) ~= sign(dx([2:end end]))];
    dxx = [0; diff(swin,2); 0];
    is_top = (swin > 0) & is_signchange & (dxx < 0);
    is_bot = (swin < 0) & is_signchange & (dxx > 0);
    spikes = is_big & (is_top | is_bot);

    % skip if insignificant spike count
    if nnz(spikes) <= param.nspikes, continue, end

    % determine median interval and classify as fast/slow ripple
    tt = i/fs:1/fs:min(i+step-1,numel(sig))/fs;
    interval = median(diff(find(spikes))) / fs;
    if param.slow(1) <= interval && interval <= param.slow(2)
      slow(end+1) = i;
    elseif param.fast(1) <= interval && interval <= param.fast(2)
      fast(end+1) = i;
    end
  end
end
