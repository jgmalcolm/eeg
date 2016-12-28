function [x_training y_training] = gather_training(data_struct, biomarker, s_window, window, offset)

% This function uses data to generate a model representing a first order
% difference equation of how stimulation affect different biomarkers. The
% main argument for this model is a directory consisting of neural
% recordings extracted from a fixed interval stimulation experiment. The
% function then calculates the neural state and resulting change in
% biomarker.


x_training              = [];
y_training              = [];

sampling_frequency      = 2000;
n_segments              = size(data_struct.model_data,1);

state_offset            = .1; % ???

for c1 = 1:size(data_struct.model_data,1)

    stimulation_duration    = data_struct.x_data(c1,1);
    stimulation_amplitude   = data_struct.x_data(c1,2);
    stimulation_frequency   = data_struct.x_data(c1,3);
    stimulation_time        = data_struct.time_data(c1,:);

    data                    = squeeze(data_struct.model_data(c1,:,:));

    % extract state segment
    state_segment           = extract_state_segment(data, window, stimulation_time, stimulation_duration, state_offset, sampling_frequency);

    % extract effect segment
    effect_segment          = extract_effect_segment(data, window, offset, stimulation_time, stimulation_duration, sampling_frequency);

    % extract_state_biomarker
    state_biomarker         = calculate_biomarker(state_segment, biomarker, s_window, sampling_frequency);

    % extract effect biomarker
    effect_biomarker        = calculate_biomarker(effect_segment, biomarker, s_window, sampling_frequency);

    y_training              = [y_training; effect_biomarker / state_biomarker];
    x_training              = [x_training;  data_struct.x_data(c1,:) ];
end

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function remove_outliers(threshold, gp_model)

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function state_segment = extract_state_segment(data, window, stimulation_time, stimulation_duration, offset, sampling_frequency)

if window == 0
   window = stimulation_duration;
end

segment_start_index = floor((stimulation_time-window-offset)*sampling_frequency);
segment_end_index   = segment_start_index + floor((window-offset)*sampling_frequency);

state_segment       = data(:, segment_start_index:segment_end_index);

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function effect_segment = extract_effect_segment(data, window, offset,  stimulation_time, stimulation_duration, sampling_frequency)

if window == 0;
    window = stimulation_duration;
end

segment_start_index = floor((stimulation_time+stimulation_duration+offset)*sampling_frequency);
segment_end_index   = segment_start_index + floor(window*sampling_frequency);

effect_segment      = data(:, segment_start_index:segment_end_index);

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function [biomarker_value] = calculate_biomarker(data, biomarker, s_window, sampling_frequency)

switch biomarker(1:3)

    case 'psd'
        [biomarker_value] = calculate_power_biomarker(data, biomarker, s_window, sampling_frequency);
    case 'coh'
        [biomarker_value] = calculate_coherence_biomarker(data, biomarker, s_window, sampling_frequency);
    case 'lin'
        [biomarker_value] = calculate_line_length_biomarker(data);
end

end

%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function biomarker_value = calculate_power_biomarker(data, biomarker, s_window, sampling_frequency)

params.Fs       = sampling_frequency;
params.tapers   = [3 4];
params.fpass    = s_window;

data            = (data - repmat(mean(data,2),1, size(data,2)));% ./ repmat(std(data,[],2), 1, size(data,2));
[S, ~]          = mtspectrumc(data',params);


if size(S,2) == 1
    biomarker_value       = mean(S);
else
    biomarker_value       = mean(sum(S));
end

end


function biomarker_value = calculate_line_length_biomarker(data)

biomarker_value = mean(sum(abs(diff(data(5:8, :)'))));

end
%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%
function biomarker = calculate_coherence_biomarker(data, biomarker, sampling_frequency)

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
% a   = a + sum(mscohere(data(1,:)',data(2,:)', [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(1,:)',data(3,:)', [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(1,:)',data(4,:)', [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(2,:)',data(3,:)', [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(2,:)',data(4,:)', [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(3,:)',data(4,:)', [], [],fpass,sampling_frequency));

a   = a + sum(mscohere(data(5,:)',data(6,:), [], [],fpass,sampling_frequency));
a   = a + sum(mscohere(data(5,:)',data(7,:), [], [],fpass,sampling_frequency));
a   = a + sum(mscohere(data(5,:)',data(8,:), [], [],fpass,sampling_frequency));
a   = a + sum(mscohere(data(6,:)',data(7,:), [], [],fpass,sampling_frequency));
a   = a + sum(mscohere(data(6,:)',data(8,:), [], [],fpass,sampling_frequency));
a   = a + sum(mscohere(data(7,:)',data(8,:), [], [],fpass,sampling_frequency));

% a   = a + sum(mscohere(data(1,:)',data(8,:), [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(2,:)',data(7,:), [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(3,:)',data(6,:), [], [],fpass,sampling_frequency));
% a   = a + sum(mscohere(data(4,:)',data(5,:), [], [],fpass,sampling_frequency));

% d   = sum(mscohere(data(4,:)',data(5,:), [], [],fpass,sampling_frequency));
q   =a/6;


biomarker = q';

end
