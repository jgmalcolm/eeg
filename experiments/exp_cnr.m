clear

KL = loadcached('divergence');
BG = loadcached('divergence_bg');

windows = unique(KL.Window);  nwindows = numel(windows);
offsets = unique(KL.Offset);  noffsets = numel(offsets);

[mu mu_bg st_bg] = deal(nan(nwindows,noffsets));
for i = 1:nwindows
  win = windows(i);
  for j = 1:noffsets
    off = offsets(j);
    ind = find(KL.Window == win & KL.Offset == off);
    div = KL(ind,:).Divergence;
    mu(i,j) = mean(div);
    mu_bg(i,j) = BG(BG.Window == win,:).Mean;
    st_bg(i,j) = BG(BG.Window == win,:).Std;
  end
end


clf
sp(1,2,1)
surf(offsets, windows, mu)
xlabel('offset')
ylabel('window')
zlabel('divergence')
set(gca,'FontSize',24)
axis vis3d

sp(1,2,2)
surf(offsets, windows, (mu-mu_bg)./st_bg)
xlabel('offset')
ylabel('window')
zlabel('CNR')
set(gca,'FontSize',24)
axis vis3d
