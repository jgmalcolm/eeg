clear

[signal stim fs] = loadcached('data/ARN045_baseline_stim_4V_2000hz', 'signal', 'stim', 'fs');
nch = size(signal,1);

%% extract pre/stim/post segments
sig_pre  = signal(:,fix(0.5*fs):fix(1.5*fs),:);
sig_stim = signal(:,fix(2.5*fs):fix(3.5*fs),:);

% sig_pre = zscore(sig_pre,0,2);
% sig_stim = zscore(sig_stim,0,2);

maxorder = 25; % maximum model order to consider
acmaxlags = 0; % estimate from data

trials = unique(stim(:,{'Channel' 'Frequency'}));
[gc_pre gc_stim] = deal([]);
for i = 1:height(trials)
  trials(i,:)
  ch = trials{i,'Channel'};
  fr = trials{i,'Frequency'};
  idx = find(stim.Channel==ch & stim.Frequency==fr);
  gc_pre(:,:,i) = est_granger(sig_pre(:,:,idx), maxorder, acmaxlags);
  gc_stim(:,:,i) = est_granger(sig_stim(:,:,idx), maxorder, acmaxlags);
end
return

% recording channels
lookup = [2 4 6 8 9 11 13 15];

channels = unique(stim.Channel);
nch = numel(channels);
freqs = unique(stim.Frequency);
nfreqs = numel(freqs);

lbl = {};
for i = 1:numel(freqs)
  lbl{i} = sprintf('%d Hz',freqs(i));
end

splookup = [1 3 5 7 8 6 4 2]; % lookup for subplots (CA3 left, CA1 right)
stim_lookup = [1 3 5 7 10 12 14 16];
read_lookup = [2 4 6 8 9 11 13 15];

rev = @(x) x(end:-1:1);

clf
for k = 1:nch
  ch = channels(k);

  clear mu se
  for i = 1:nfreqs
    freq = freqs(i);

    idx = find(stim.Channel==ch & stim.Frequency==freq);
    change = Sstim(:,:,idx) - Spre(:,:,idx);
    mu(:,:,i) = mean(change,3);
    se(:,:,i) = 2*stderr(change,3);
  end

  ff_ = [ff rev(ff)];
  for i = 1:nch
    sp(nch,nch,(k-1)*nch + i)
    mu_ = squeeze(mu(:,i,:));
    se_ = squeeze(se(:,i,:));
    clear h
    cc = [1 0 0; 0 .7 0; 0 0 1];
    for j = 1:nfreqs
      m = mu_(:,j)';
      s = se_(:,j)';
      h(j) = patch(ff_, ...
                   [m+s rev(m-s)], ...
                   cc(j,:), ...
                   'EdgeColor','none', 'FaceAlpha',.3);
    end
    title(sprintf('stim %d, read %d', stim_lookup(k), read_lookup(i)))
    axis tight
    set(gca, 'FontSize',12)
    if k==1 && i==1
      legend(h, lbl{:}, 'Location','Best')
    end
  end

end

% use same YLim for all subplots
children = get(gcf, 'Children');
lim = [];
for i = 1:numel(children)
  if isa(children(i),'matlab.graphics.axis.Axes')
    lim(:,end+1) = get(children(i),'YLim');
  end
end
lim = [min(lim(:)) max(lim(:))];
for i = 1:numel(children)
  if isa(children(i),'matlab.graphics.axis.Axes')
    set(children(i),'YLim',lim);
  end
end
