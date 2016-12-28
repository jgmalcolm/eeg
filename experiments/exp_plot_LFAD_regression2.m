function dloc = exp_plot_LFAD_regression2(data, model, model_args)

  kfolds = 30;
  rand('state',0) % determinism

  locations = unique(data.Location)';
  nlocs = numel(locations);

  dloc = table();
  for loc = 1:nlocs
    loc = locations(loc);

    % gather data
    d = data(data.Location==loc & data.Amplitude > 0 & data.Frequency > 0,:);
    X = d(:,{'Biomarker' 'Frequency' 'Amplitude' 'Duration' 'dxBiomarker'});
    Y = d.dxBiomarker(:);

    % determine folds used for all models
    n = numel(Y);
    cv = cvpartition(n,'KFold',kfolds);

    preds = {'Frequency' 'Amplitude' 'Duration'};

    %% compare models
    args = {'ResponseVar','dxBiomarker', 'PredictorVars'};
    bio_mdl     = @(X,x) predict(feval(model, X, model_args{:}, args{:}, {'Biomarker'}), x);
    stim_mdl    = @(X,x) predict(feval(model, X, model_args{:}, args{:}, preds), x);
    biostim_mdl = @(X,x) predict(feval(model, X, model_args{:}, args{:}, {'Biomarker' preds{:}}), x);
    [Y_bio mse_bio] = loocv(X,Y,bio_mdl);
    [Y_stim  mse_stim] = loocv(X,Y,stim_mdl);
    [Y_biostim mse_biostim] = loocv(X,Y,biostim_mdl);
    fprintf('stim %.5f   bio %.5f   biostim %.5f    %-s\n', mse_stim, mse_bio, mse_biostim, char(loc))

    dloc = [dloc; table(loc, n, {Y}, {Y_stim}, {Y_bio}, {Y_biostim})];
  end

  vars = {'Location' 'N' 'Y' 'Y_stim' 'Y_bio' 'Y_biostim'};
  dloc.Properties.VariableNames = vars;
  return

  % final table sorted by MSE of biostim
  mse_bio = cellfun(@(x,y) mean((x-y).^2), dloc.Y, dloc.Y_bio);
  mse_stim = cellfun(@(x,y) mean((x-y).^2), dloc.Y, dloc.Y_stim);
  mse_biostim = cellfun(@(x,y) mean((x-y).^2), dloc.Y, dloc.Y_biostim);
  nn = cellfun(@numel, dloc.Y);

  dloc = table(dloc.Location, nn, mse_stim, mse_bio, mse_biostim);
  dloc.Properties.VariableNames = {'Location' 'N' 'mseStim' 'mseBio' 'mseBioStim'};
  [~,idx] = sort(mse_biostim); % sort by biostim
  dloc = dloc(idx,:); % sort by biostim
end


% k-fold cross-validation
function [y_hat mse] = kfoldcv(x,y,predict,cv)
  n = size(x,1);

  y_hat = nan(n,1);
  for i = 1:cv.NumTestSets
    test = find(cv.test(i));
    train = find(cv.training(i));
    x_test = x(test,:);
    x_train = x(train,:);

    y_hat(test) = predict(x_train, x_test);
  end
  mse = mean((y - y_hat).^2);
end



% leave one out cross-validation
function [y_hat mse] = loocv(x,y,predict)
  n = size(x,1);
  y_hat = nan(n,1);
  idx = 1:n;
  parfor i = idx
    y_hat(i) = predict(x(idx~=i,:), x(i,:));
  end
  mse = mean((y - y_hat).^2);
end
