function MI = est_mutualinfo(signal)
% Estimate pairwise channel mutual information
%
%  SIGNAL channels x time
%  MI pairwise mutual information (diagonal=1)

  if ~exist('mutualinfo')
    addpath('~/src/mutualinfo')
  end

  nch = size(signal,1);
  signal_ = zscore(signal,0,2);
  MI = zeros(nch);
  for i = 1:nch
    si = signal_(i,:);
    for j = i+1:nch
      sj = signal_(j,:);
      [MI(i,j) MI(j,i)] = deal(mutualinfo(si,sj));
    end
  end
  MI(1:size(MI,1)+1:end) = 1; %% diagonal self==1
end
