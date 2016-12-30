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





function [state_segment idx] = extract_state_segment(data, window, stimulation_time, stimulation_duration, offset, sampling_frequency)

  segment_start_index = floor((stimulation_time-window-offset)*sampling_frequency);
  segment_end_index   = segment_start_index + floor(window*sampling_frequency);

  idx = segment_start_index:segment_end_index;
  state_segment       = data(:, idx);

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function [effect_segment idx] = extract_effect_segment(data, window, offset,  stimulation_time, stimulation_duration, sampling_frequency)

segment_start_index = floor((stimulation_time+stimulation_duration+offset)*sampling_frequency);
segment_end_index   = segment_start_index + floor(window*sampling_frequency);

idx = segment_start_index:segment_end_index;
effect_segment      = data(:, idx);

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function [biomarker_value] = calculate_biomarker(data, biomarker, freqs, ds)

switch biomarker(1:3)

    case 'psd'
        [biomarker_value] = calculate_power_biomarker(data, biomarker, freqs, ds);
    case 'coh'
        [biomarker_value] = calculate_coherence_biomarker(data, biomarker, freqs, ds);
    case 'lin'
        [biomarker_value] = calculate_line_length_biomarker(data);
end

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function biomarker_value = calculate_power_biomarker(data, biomarker, freqs, ds)

params.Fs       = ds;
params.tapers   = [3 4];
params.fpass    = freqs;

% data            = (data - repmat(mean(data,2),1, size(data,2)));% ./ repmat(std(data,[],2), 1, size(data,2)); % ???
data = zscore(data,0,2);
[S, ~] = mtspectrumc(data',params);


if size(S,2) == 1
    biomarker_value       = mean(S);
else
    biomarker_value       = mean(sum(S));
end

end


function biomarker_value = calculate_line_length_biomarker(data)

biomarker_value = mean(sum(abs(diff(data(5:8, :)'))));

end


function biomarker = calculate_xcorr_biomarker(data)
% @(x,y) max(abs(xcorr(x,y))) / (norm(x,2)*norm(y,2))
end


%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function biomarker = calculate_coherence_biomarker(data, biomarker, ds)

switch biomarker

    case 'coh_theta'
        fpass = 4:.1:7;
    case 'coh_gamma'
        fpass = 32:.1:50;
    case 'coh_alpha'
        fpass = 8:.1:15;
    case 'coh_beta'
        fpass = 16:.1:31;
    case 'coh_all'
        fpass = 1:.1:50;

end

a   = 0;
% a   = a + sum(mscohere(data(1,:)',data(2,:)', [], [],fpass,ds));
% a   = a + sum(mscohere(data(1,:)',data(3,:)', [], [],fpass,ds));
% a   = a + sum(mscohere(data(1,:)',data(4,:)', [], [],fpass,ds));
% a   = a + sum(mscohere(data(2,:)',data(3,:)', [], [],fpass,ds));
% a   = a + sum(mscohere(data(2,:)',data(4,:)', [], [],fpass,ds));
% a   = a + sum(mscohere(data(3,:)',data(4,:)', [], [],fpass,ds));

a   = a + sum(mscohere(data(5,:)',data(6,:), [], [],fpass,ds));
a   = a + sum(mscohere(data(5,:)',data(7,:), [], [],fpass,ds));
a   = a + sum(mscohere(data(5,:)',data(8,:), [], [],fpass,ds));
a   = a + sum(mscohere(data(6,:)',data(7,:), [], [],fpass,ds));
a   = a + sum(mscohere(data(6,:)',data(8,:), [], [],fpass,ds));
a   = a + sum(mscohere(data(7,:)',data(8,:), [], [],fpass,ds));

% a   = a + sum(mscohere(data(1,:)',data(8,:), [], [],fpass,ds));
% a   = a + sum(mscohere(data(2,:)',data(7,:), [], [],fpass,ds));
% a   = a + sum(mscohere(data(3,:)',data(6,:), [], [],fpass,ds));
% a   = a + sum(mscohere(data(4,:)',data(5,:), [], [],fpass,ds));

% d   = sum(mscohere(data(4,:)',data(5,:), [], [],fpass,ds));
q   =a/6;


biomarker = q';

end
