clear
% biomarker_data
data = loadcached('data/stim_par_table','data');

subjects = unique(data.Subject);
locations = unique(data.Location);
preds = {'Frequency' 'Amplitude' 'Duration'};

% last column is across all subjects
[rr2 pvals] = deal(nan(numel(locations), numel(subjects)+1));
%% for each subject
for subj = unique(data.Subject)'
  data_subj = data(data.Subject==subj,:);

  %% for each location of subject
  locs = unique(data_subj.Location);
  for loc = 1:numel(locs)
    d = data_subj(data_subj.Location==locs(loc),:);

    % use predictors with more than one unique value (ignore NaN)
    preds_ = preds(cellfun(@(p) numel(unique(notnan(d.(p)))) > 1, preds));

    lm = fitlm(d,'linear', ...
               'PredictorVars',{preds_{:} 'Biomarker'}, ...
               'ResponseVar','dxBiomarker');
    idx = locations == locs(loc);
    rr2(idx,subj) = lm.Rsquared.Adjusted;
    pvals(idx,subj) = coefTest(lm);
  end
end


% last column is across all subjects
for loc = 1:numel(locations)
  d = data(data.Location==locations(loc),:);

  % use predictors with more than one unique value (ignore NaN)
  preds_ = preds(cellfun(@(p) numel(unique(notnan(d.(p)))) > 1, preds));

  lm = fitlm(d,'linear', ...
             'PredictorVars',{preds_{:} 'Biomarker'}, ...
             'ResponseVar','dxBiomarker');
  rr2(loc,end) = lm.Rsquared.Adjusted;
  pvals(loc,end) = coefTest(lm);
end

% reorder by r^2
[~,idx] = sort(rr2(:,end),1,'descend');
rr2_ = rr2(idx,:);
pvals_ = pvals(idx,:);
locations_ = locations(idx);

colormap jet
h = imagesc(rr2_);
set(h, 'AlphaData', ~isnan(rr2_)) % white NaN
grid on
xlabel('subject')
colorbar('WestOutside')
set(gca, 'FontSize',14, ...
         'YTick', 1:numel(locations), ...
         'YTickLabel', char(locations_), ...
         'YAxisLocation', 'right')
xx = get(gca, 'XTickLabel');
xx{end} = 'all';
set(gca, 'XTickLabel', xx)

%{
% indicate significant locations
% [pvals_ idx] = bonf_holm(pvals, 0.05);
[ll ss] = find(pvals_ < 0.05);
hold on
plot(ss,ll,'w*', 'MarkerSize',12, 'LineWidth',2)
hold off
%}



return
for i = 1:size(X,2)
  sp(2,2,i)
  plot(X(:,i), Y, '.')
  grid on, axis tight
end

%{
for i = 1:3
  x = X(:,i+1);
  xx = unique(x);
  [y_mu y_sd] = deal(nan(size(xx)));
  for j = 1:numel(xx)
    y = Y(x == xx(j));
    y_mu(j) = mean(y);
    y_sd(j) = stderr(y);
  end
  sp(nloc,3,(loc-1)*3+i)
  errorbar(xx,y_mu,y_sd, 'LineWidth',3,'Color','b');
  grid on, axis tight, box off
end
%}

%{
nsamples = numel(Y);
Y_linear = cell(nloc,1);
for loc = 1:nloc
fprintf('loc %d\n', loc)

clear y_*
parfor i = 1:nsamples
  % remove this sample from training set
  [X_ Y_] = deal(X,Y);
  X_(i,:) = [];
  Y_(i) = [];
  x = X(i,:);

  %{
  % lasso
  [B fitinfo] = lasso(X_,Y_, 'CV',10);
  idx = fitinfo.Index1SE;
  b = B(:,idx);
  b0 = fitinfo.Intercept(idx);
  y_lasso(i) = x * b + b0;
  %}

  %{
  % neural net
  nhidden = 20; % can be high because using regularization
  nnet = feedforwardnet(nhidden,'trainbr');
  nnet.trainParam.showWindow = false;
  nnet = train(nnet, X_', Y_');
  y_nn(i) = nnet(x');
  %}

  %{
  % regression tree
  tree = fitrtree(X_, Y_);
  y_tree(i) = predict(tree, x);
  %}

  %{
  % svm - gaussian
  svm = fitrsvm(X_,Y_, 'KernelFunction','gaussian','KernelScale','auto','Standardize',true);
  y_svm_gauss(i) = predict(svm, x);
  %}

  % linear
  b = [X_ ones(nsamples-1,1)] \ Y_;
  y_linear(i) = [x 1] * b;
end


Y_linear{loc} = y_linear;
end

[r2 rho] = deal(nan(size(bio_diff)));
stim_Y = [];
for loc = 1:nloc

  bio = bio_pre{loc}(:);
  stim = double(stim_param{loc});
  change = bio_diff{loc}(:);

  % if loc <= 26
  %   stim = stim(:,[1 3]); % frequency, duration
  % elseif loc > 26
  %   stim = stim(:,[1 2]); % frequency, amplitude
  % end
  % remove -999 stimulations
  [rows cols] = find(stim == -999);
  Y = change;
  Y(rows) = [];
  stim(rows,:) = [];

  Y_hat = Y_linear{loc}(:);

  stim_Y = [stim_Y; stim Y Y_hat];
end

%% freq
x = stim_Y(:,1);
show_r2_rho(x,stim_Y(:,4),stim_Y(:,5),'Hz')

%% amplitude
x = stim_Y(:,2);
show_r2_rho(x,stim_Y(:,4),stim_Y(:,5),'mV')

%% duration
x = stim_Y(:,3);
show_r2_rho(x,stim_Y(:,4),stim_Y(:,5),'ms')

return

y = stim_se(:,4);
xx = unique(x);
fprintf('\n%5s  %10s  %10s\n', 'mV', 'median', 'iqr')
for i = 1:numel(xx)
  yy = y(x == xx(i));
  fprintf('%5g  %10.5f  %10.5f\n', ...
          xx(i), median(yy), iqr(yy))
end

%% duration
x = stim_se(:,3);
y = stim_se(:,4);
xx = unique(x);
fprintf('\n%5s  %10s  %10s\n', 'ms', 'median', 'iqr')
for i = 1:numel(xx)
  yy = y(x == xx(i));
  fprintf('%5g  %10.5f  %10.5f\n', ...
          xx(i), median(yy), iqr(yy))
end
%}
