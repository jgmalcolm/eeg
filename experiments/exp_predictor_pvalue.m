function exp_predictor_pvalue(data, locs)

  kfolds = 30;
  rand('state',0) % determinism

  preds = {'Biomarker' 'Frequency' 'Amplitude' 'Duration'};
  args = {'ResponseVar','dxBiomarker', 'PredictorVars', preds};

  hold on
  h = [];
  off = [-.05 0 .05];
  for i = 1:numel(locs)
    % gather data
    d = data(data.Location==locs{i} & data.Amplitude > 0 & data.Frequency > 0,:);
    X = d(:,{'Biomarker' 'Frequency' 'Amplitude' 'Duration' 'dxBiomarker'});
    Y = d.dxBiomarker(:);

    % determine folds to use
    n = numel(Y);
    cv = cvpartition(n,'KFold',kfolds);

    % gather and plot p-values
    fitmodel = @(x_train) fitlm(x_train, 'linear', args{:});
    pvalues = kfoldcv(X,cv,fitmodel);
    h(i) = errorbar((1:4)+off(i), mean(pvalues,2), std(pvalues,0,2)/sqrt(kfolds));
  end

  plot([0 4]+.5, [0.05 0.05],'k--');
  hold off
  set(h, 'LineStyle','none', 'LineWidth',3, 'MarkerSize',20, 'Marker','.');
  set(gca, 'FontSize',20, ...
           'XTick', 1:4, 'XTickLabel', preds, ...
           'YTick', [0.05 .1:.1:max(ylim)])
  ylabel('p-value')
  title('p-value of model terms')

  legend(h, locs{:}, 'Location','Best')
end


% k-fold cross-validation
function pvalue = kfoldcv(x,cv,fitmodel)
  n = size(x,1);

  pvalue = [];
  for i = 1:cv.NumTestSets
    train = find(cv.training(i));
    mdl = fitmodel(x(train,:));
    pvalue = [pvalue mdl.Coefficients.pValue(2:end)];
  end
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
