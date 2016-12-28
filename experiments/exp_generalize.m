function exp_generalize(data)

  rand('state',0) % determinism


  [pos_lin pos_qua pos_svm] = calc_mse(data, 'Sub');
  [neg_lin neg_qua neg_svm] = calc_mse(data, 'TC');

  % [pos_lin pos_qua pos_svm] = deal(.7, .8, .9);
  % [neg_lin neg_qua neg_svm] = deal(.8, .9, 1.);

  X = [.1 .2 .3;  .7 .8 .9];
  Y = [pos_lin pos_qua pos_svm;
       neg_lin neg_qua neg_svm];

  h = bar(X,Y);
  set(h, 'BarWidth', 2.5, 'EdgeAlpha',0)
  set(h(1), 'FaceColor', hex2rgb('418CF0')/256)
  set(h(2), 'FaceColor', hex2rgb('FCB441')/256)
  set(h(3), 'FaceColor', hex2rgb('2E8B57')/256)

  set(gca, 'FontSize',20, ...
           'YGrid', 'on', ...
           'XTick', [.2 .8], 'XTickLabel', {'Sub' 'TC'}, 'XLim', [0 1])
  ylabel('MSE')
  title('Generalizability across models')
  legend(h, 'Linear', 'Quadratic', 'SVM', 'Location','South')
end


% leave one out cross-validation
function mse = loocv(x,y,predict)
  n = size(x,1);
  y_hat = nan(n,1);
  idx = 1:n;
  parfor i = idx
    y_hat(i) = predict(x(idx~=i,:), x(i,:));
  end
  mse = mean((y - y_hat).^2);
end


function [lin qua svm] = calc_mse(data, loc)

  preds = {'Biomarker' 'Frequency' 'Amplitude' 'Duration'};
  args = {'ResponseVar','dxBiomarker', 'PredictorVars', preds};

  % gather data
  d = data(data.Location==loc & data.Amplitude > 0 & data.Frequency > 0,:);
  X = d(:,{'Biomarker' 'Frequency' 'Amplitude' 'Duration' 'dxBiomarker'});
  Y = d.dxBiomarker(:);

  % warning('using fraction of data for testing')
  % X = X(1:100,:);
  % Y = Y(1:100);

  % form the models
  lin = @(X,x) predict(fitlm(X,'linear',args{:}), x);
  qua = @(X,x) predict(fitlm(X,'quadratic',args{:}), x);
  svm = @(X,x) predict(fitrsvm(X,'dxBiomarker', ...
                               'PredictorNames', preds, ...
                               'ResponseName','dxBiomarker', ...
                               'Standardize',true, ...
                               'KernelFunction','gaussian', ...
                               'KernelScale','auto'), x);

  % calculate MSE using leave-one-out
  lin = loocv(X,Y,lin)
  qua = loocv(X,Y,qua)
  svm = loocv(X,Y,svm)
end
