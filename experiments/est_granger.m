function [GC moaic ac] = est_granger(signal, maxorder, acmaxlags)
% Estimate the pairwise conditional Granger Causality
%
% [GC MO AC] = EST_GRANGER(SIGNAL, MAXORDER, ACMAXLAGS)
%
% Input:
%   SIGNAL [channels x observations]
%   MAXORDER maximum model order to consider
%   ACMAXLAGS maximum autocorr lags to use (0 = calculate from data)
% Output:
%   GC - pairwise conditional Granger Causality (zero diagonal)
%   MO - model order by Akaike info criterion
%   AC - autocorrelation lags

  if ~exist('tsdata_to_infocrit')
    run '~/src/mvgc_v1.0/startup.m'
  end

  ac = nan;
  nch = size(signal,1);
  [aic bic moaic mobic] = tsdata_to_infocrit(signal, maxorder, 'LWR', false);
  [A SIG] = tsdata_to_var(signal, moaic); % use AIC
  if isbad(A)
    warning('VAR estimation failed')
    GC = nan(nch);
    return
  end
  [AC info] = var_to_autocov(A, SIG, acmaxlags);
  if info.error, warning(info.errmsg), end
  for i=1:info.warnings, warning(info.warnmsg{i}), end
  GC = autocov_to_pwcgc(AC);
  if isbad(GC,false)
    warning('GC calculation failed')
    GC = nan(nch);
    return
  end
  GC(1:size(GC,1)+1:end) = 0; %% diagonal self == 0
  ac = info.aclags;
end
