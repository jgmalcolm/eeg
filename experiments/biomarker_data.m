clear
[bio_pre stim_param bio_diff locs] = loadcached('data/stim_par','Prob_pre','Stim_par','Prob_dif','Stim_loc');

nsubj = numel(bio_pre);

data = table();
for subj = 1:nsubj
  bio = bio_pre{subj}(:);
  stim = double(stim_param{subj});
  stim(stim == -999) = 0; % -999 duration ==> Pulse
  change = bio_diff{subj}(:);

  data_ = table(subj*ones(size(bio)), bio, stim(:,1), stim(:,2), stim(:,3), locs{subj}', change);
  data = [data; data_]; % append rows
end

vars = {'Subject' 'Biomarker' 'Frequency' 'Amplitude' 'Duration' 'Location' 'dxBiomarker'};
data.Properties.VariableNames = vars;
data.Subject = categorical(data.Subject);
data.Location = categorical(data.Location);

save('data/stim_par_table','data')
