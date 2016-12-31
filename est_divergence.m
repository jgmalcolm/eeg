function div = est_divergence(X,Y)
% estimate of KL divergence via nearest-neighbor method (Wang2006)
% X,Y - each row is a channel

  % force to double
  X = double(X);
  Y = double(Y);

  % orient to row vectors
  if iscolumn(X), X = X'; end
  if iscolumn(Y), Y = Y'; end

  % check dimensions
  assert(ndims(X) == 2)
  assert(ndims(Y) == 2)
  [nch nx] = size(X);
  assert(nch == size(Y,1)) % same channel count
  ny = size(Y,2);

  % acc = 0;
  % pts = 0;
  % for i = 1:nx
  %   x = X(:,i);
  %   dxx2 = sort(sum(bsxfun(@minus,X,x).^2,1));
  %   dxx2 = dxx2(2); % skip self
  %   if dxx2 > 0 && ~isinf(dxx2)
  %     dxy2 = min(sum(bsxfun(@minus,Y,x).^2,1));
  %     if dxy2 > 0 && ~isinf(dxy2)
  %       acc = acc + log(dxy2 / dxx2);
  %       pts = pts + 1;
  %     end
  %   end
  % end
  % acc = acc / pts;
  acc = est_divergencemex(X,Y);

  div = (1/2)*nch*acc + log(ny/(nx-1));

  % (1/2) is pulled off norm inside loop via log
end


% "A Nearest-Neighbor Approach to Estimating Divergence between
% Continuous Random Vectors." Qing Wang, Sanjeev R. Kulkarni, Sergio
% Verdu, IEEE International Symposium on Information Theory (ISIT),
% 2006.
