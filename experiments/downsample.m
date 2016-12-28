function [s_ fs_] = downsample(s,fs,fs_)
% Downsample a signal
% [S_ FS_] = DOWNSAMPLE(S,FS,FS_)
%  S,FS input signal (channels,time) and sample frequency
%  S_,FS_ output signal and sample frequency

  assert(fs_ < fs, 'output sampling rate (%d Hz) must be lower than input (%d Hz)', fs_,fs)

  [nch dt ntrials] = size(s);
  % permute time to first dimension
  s = permute(s, [2 1 3]);

  % do first vector to see result size
  s_ = resample(double(s(:,1)), fs_, fs);
  s_(:,nch,ntrials) = 0; % preallocate

  for i = 2:nch*ntrials
    s_(:,i) = resample(double(s(:,i)), fs_, fs);
  end

  s_ = permute(s_, [2 1 3]);
