function progress = parforprogress(N_)
  persistent fname N

  if nargin
    N = N_;
    id = feature('getpid'); % process ID of matlab
    fname = sprintf('/tmp/parfor-%d', id);
    fid = fopen(fname, 'w');
    % fprintf(fid, '%d\n', N)
    % fwrite(fid, N,'integer*4')
    assert(fclose(fid) == 0);

  else
    fid = fopen(fname,'a+');
    assert(fseek(fid,0,-1) == 0) % reset to beginning
    status = fscanf(fid, '%d'); % read completed
    assert(fprintf(fid, '1\n') == 2) % add one

    assert(fclose(fid) == 0);

    progress = sum(status) / N;
  end
end
