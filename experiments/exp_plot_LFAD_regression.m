function exp_plot_LFAD_regression(data,loc,freq,dur,amp)

  rand('state',0) % determinism for SVM calculations

  % gather data
  d = data(data.Location==loc & data.Frequency==freq & data.Duration==dur & data.Amplitude==amp,:);
  [X idx] = sort(d.Biomarker(:));
  Y = d.dxBiomarker(idx);

  cla
  hold on
  h = plot(X,Y,'ko','MarkerSize',10);

  %% SVM
  svm_fn = @(X,Y,x) predict(fitrsvm(X,Y', ...
                                    'BoxConstraint', 0.3, ...
                                    'KernelFunction','gaussian', ...
                                    'KernelScale','auto', ...
                                    'Standardize',true), ...
                            x');
  [Y_svm mse_svm] = loocv(X,Y,svm_fn);
  h(end+1) = plot(X,Y_svm,'r.','MarkerSize',20);

  %% Linear
  lm_fn = @(X,Y,x) predict(fitlm(X,Y','linear'), x');
  [Y_lm mse_lm] = loocv(X,Y,lm_fn);
  h(end+1) = plot(X,Y_lm,'b.','MarkerSize',20);

  %% Polynomial
  poly_fn = @(X,Y,x) predict(fitlm(X,Y','poly3'), x');
  [Y_poly mse_poly] = loocv(X,Y,poly_fn);
  h(end+1) = plot(X,Y_poly,'m.','MarkerSize',20);

  legend(h, ...
         'observations', ...
         sprintf('svm regression (mse %.4f)', mse_svm), ...
         sprintf('linear regression (mse %.4f)', mse_lm), ...
         sprintf('polynomial regression (mse %.4f)', mse_poly), ...
         'Location','Best')

  axis tight
  xlabel('Biomarker')
  ylabel('Change in Biomarker')
  title(sprintf('%s  %d Hz  %d ms  %.2f mA', loc, freq, dur, amp))

  set(gca, 'FontSize',16, 'LineWidth',1)
  hold off
end

function y_hat = pls_fn(X,Y,x)
  [~,~,~,~,beta] = plsregress(X,Y,1);
  y_hat = [ones(size(x,1),1) x] * beta;
end

function [bc fval] = opt_hyperparam(X,Y)
  svm_fn = @(bc) loocv(X,Y, ...
                       @(X,Y,x) predict(fitrsvm(X,Y', ...
                                                'BoxConstraint', bc, ...
                                                'KernelFunction','gaussian', ...
                                                'KernelScale','auto', ...
                                                'Standardize',true), ...
                                        x'));
  x0 = 0.3;
  opts = optimoptions(@fmincon, 'Display','iter');
  [bc fval] = fmincon(svm_fn, x0, [],[],[],[], 0,5, [], opts)
end

function [y_hat mse] = loocv(xx,yy, mdl_fn)
  y_hat = nan(size(xx));
  for i = 1:numel(xx)
    % leave one out
    [xx_ yy_] = deal(xx,yy);
    xx_(i) = [];
    yy_(i) = [];

    y_hat(i) = mdl_fn(xx_, yy_, xx(i));
  end
  mse = mean((yy - y_hat).^2);
end
