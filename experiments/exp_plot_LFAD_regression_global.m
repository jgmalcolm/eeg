function dloc = exp_plot_LFAD_regression_global(data)

  kfolds = 30;
  rand('state',0) % determinism

  dloc = {};
  locations = unique(data.Location)';
  nlocs = numel(locations);

  % compare against each location
  parfor i = 1:nlocs
    loc = locations(i);

    % gather global data from other locations
    X_other = data(data.Location~=loc & data.Amplitude > 0 & data.Frequency > 0,:);

    % gather local data
    X = data(data.Location==loc & data.Amplitude > 0 & data.Frequency > 0,:);
    Y = X.dxBiomarker(:);

    % determine predictor variables
    preds = {'Frequency' 'Amplitude' 'Duration'};
    preds_ = preds(cellfun(@(p) numel(unique(X.(p))) > 1, preds));
    assert(numel(preds_) > 0) % intend to regress against more than just biomarker

    % determine folds
    n = numel(Y);
    cv = cvpartition(n,'KFold',kfolds);

    % compare local vs global
    args = {'ResponseVar','dxBiomarker', 'PredictorVars'};
    local_mdl  = @(X,x) predict(fitlm(X, 'quadratic', args{:}, {'Biomarker' preds{:}}), x);
    global_mdl = @(X,x) predict(fitlm([X;X_other], 'quadratic', args{:}, {'Biomarker' preds{:}}), x);
    [~,mse_local]  = kfoldcv(X,Y,local_mdl,cv);
    [~,mse_global] = kfoldcv(X,Y,global_mdl,cv);
    fprintf('local %.5f   global %.5f     %-s\n', mse_local, mse_global, char(loc))

    dloc(i,:) = {loc, n, mse_local, mse_global};
  end

  vars = {'Location' 'N' 'mseLocal' 'mseGlobal'};
  dloc = cell2table(dloc, 'VariableNames', vars);
  dloc = sortrows(dloc,'mseGlobal');
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
