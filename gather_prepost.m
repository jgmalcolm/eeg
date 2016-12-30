function [Spre Spost] = gather_prepost(data_struct, freqs, window, offset)

% data_struct
%   - x_data = [dur amp freq]
%   - time_data = stimulation times
%   - model_data = [stims * channels * time]
% stim_time: time points of stimulation (seconds)
% stim_params: [dur amp freq]
% signal: multi-channel signal (channels * time)
% freqs: frequency range [min max]
% window: how many seconds to consider as window
% offset: how soon after end of stimulation to grab window

  ds = 2000; % Hz
  state_offset = .1; % time between end of pre-stim window and start of stim (seconds)

  nsamples = size(data_struct.model_data,1);
  [Spre Spost] = deal(zeros(0,0,0));

  for i = 1:nsamples
    stim_duration  = data_struct.x_data(i,1);
    stim_amplitude = data_struct.x_data(i,2);
    stim_frequency = data_struct.x_data(i,3);
    stim_time      = data_struct.time_data(i,:);

    % extract segments
    signal = squeeze(data_struct.model_data(i,:,:));
    [spre preidx] = extract_state_segment(signal, window, stim_time, stim_duration, state_offset, ds);
    [spost postidx] = extract_effect_segment(signal, window, offset, stim_time, stim_duration, ds);
    assert(isequal(size(spre),size(spost)))

    Spre(:,:,i) = spre;
    Spost(:,:,i) = spost;
  end

end



function [segment idx] = extract_state_segment(data, window, stimulation_time, stimulation_duration, offset, sampling_frequency)

  start_idx = floor((stimulation_time-window-offset)*sampling_frequency);
  end_idx   = start_idx + floor(window*sampling_frequency);

  idx = start_idx:end_idx;
  segment = data(:, idx);

end

function [segment idx] = extract_effect_segment(data, window, offset,  stimulation_time, stimulation_duration, sampling_frequency)

  start_idx = floor((stimulation_time+stimulation_duration+offset)*sampling_frequency);
  end_idx   = start_idx + floor(window*sampling_frequency);

  idx = start_idx:end_idx;
  segment = data(:, idx);

end
