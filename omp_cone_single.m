function [x, Da, support] = omp_cone_single(y, D, s, rad, iter_cd, err_thresh)

% Optimal Matching Pursuit with cone atoms
% 
% Input:
%   y         - signal
%   D         - dictionary (normalized atoms)
%   s         - sparsity level
%   rad       - cone radii (if scalar, then all atoms have the same radius)
%   iter_cd   - number of iterations in the coordinate descent refinement
%               (default = 1)
%   err_thresh- representation error for which OMP is stopped
%
% Output:
%   X         - sparse representations matrix
%   Da        - actual dictionary for last representation (meaningful when N=1)
%   support   - indices of atoms used in representation

% BD 14.06.2022

if nargin < 5
  iter_cd = 1;
end

if nargin < 6
  err_thresh = 1e-14;
end

% sizes and constants
[m,n] = size(D);
if length(rad) == 1
  rad = rad * ones(n,1);
end
rad2 = rad.*rad/2;
lam = 1 - rad2;

x = zeros(s,1);
Da = zeros(m,s);

r = y;         % init residual
support = [];
for k = 1 : s
  rn = r / norm(r); % normalized residual, to be used in convex combination with best atom
  % find best current atom
  proj = D'*rn;         % projections, like in OMP
  proj(support) = 0;    % avoid atom duplication
  [~,jmax] = max(lam.*abs(proj) + sqrt((1-lam.*lam).*(1-proj.*proj))); % find nearest cone
  if proj(jmax) < 0
    rn = -rn;
  end
  p = abs(proj(jmax));
  if p + rad2(jmax) >= 1 % residual is inside cone
    dc = rn;             % hence new atom is the residual
    res_inside_cone = 1;
  else                   % residual outside cone
    b = sqrt(1/(1-p*p));
    q = -b*p*D(:,jmax) + b*rn;    % vector orthogonal on atom, in the plane atom-residual
    a = 1 - rad2(jmax);
    dc = a*D(:,jmax) + sqrt(1-a*a)*q;   % the new atom
    dc = dc/norm(dc);  % normalize, to be sure (although it should have norm 1)
    res_inside_cone = 0;
  end
  support = [support jmax]; % update support with new atom
    
  % compute (nearly) optimal representation
  Da(:,k) = dc;    % current dictionary
  if res_inside_cone    % with new atom, residual becomes zero
    x(k) = dc'*r;
    x = x(1:k);         % cut representation to useful part
    Da = Da(:,1:k);
    return
  else
    Dk = Da(:,1:k); 
    xk = Dk \ y;  % least square solution, with current atoms
    if k > 1
      for icd = 1 : iter_cd
        r = y - Dk*xk;
        %norm(r)/norm(y)
        for i = 1 : k       % a round of coordinate descent
          r = r + Dk(:,i)*xk(i);
          rn = r / norm(r);
          p = D(:,support(i))'*rn;
          if p<0
            rn = -rn;
          end
          p = abs(p);
          if p + rad2(support(i)) >= 1    % residual in the cone => ready
            Dk(:,i) = rn;
            xk(i) = rn'*r;
            x = xk;
            Da = Dk;
            return
          else
            b = sqrt(1/(1-p*p));
            q = -b*p*D(:,support(i)) + b*rn;
            a = 1 - rad2(support(i));
            dc = a*D(:,support(i)) + sqrt(1-a*a)*q;   % the new atom
            dc = dc/norm(dc);  % ???
            Dk(:,i) = dc;
            xk(i) = dc'*r;
            r = r - dc*xk(i);
          end
        end
      end
    end
    r = y - Dk*xk;
    err = norm(r)/norm(y);
    x(1:k) = xk;
    Da(:,1:k) = Dk;
    if err < err_thresh
      x = x(1:k);
      Da = Da(:,1:k);
      return
    end
  end
end
