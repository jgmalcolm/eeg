function dloc_nmse = prepare_nmse(dloc_y)

  nlocs = size(dloc_y,1);
  [nmse_stim nmse_bio nmse_biostim] = deal(nan(nlocs,1));
  for i = 1:nlocs
    nmse_stim(i)    = calc_nmse(dloc_y.Y{i}, dloc_y.Y_stim{i});
    nmse_bio(i)     = calc_nmse(dloc_y.Y{i}, dloc_y.Y_bio{i});
    nmse_biostim(i) = calc_nmse(dloc_y.Y{i}, dloc_y.Y_biostim{i});
  end

  dloc_nmse = table(dloc_y.Location, dloc_y.N, nmse_stim, nmse_bio, nmse_biostim, ...
                    'VariableNames',{'Location' 'N' 'nmseStim' 'nmseBio' 'nmseBioStim'});
end

function e = calc_nmse(y,y_hat)
  e = mean((y - y_hat).^2) / mean(y) / mean(y_hat);
  e = mean((y - y_hat).^2) / (mean(y) * mean(y_hat));
end
