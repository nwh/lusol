function [L U P Q] = lusol(A,pm_opt,varargin)
  % lusol  perform LU factorization with LUSOL

  % handle optional output
  if nargout == 2
    perm_flag = false;
  elseif nargout == 4
    perm_flag = true;
  else
    error('lusol:output','incorrect number of output arguments')
  end

  % handle optional input
  pm_string = 'matrix';
  if nargin >= 2 && ~isempty(pm_opt) && strcmpi(pm_opt,'vector')
    pm_string = 'vector';
  end

  % perform the factorization
  lu_obj = lusol_obj(A,varargin{:});

  % extract matrix factors
  if perm_flag
    % permutation vectors/matrices are requested
    [L p] = lu_obj.L0(pm_string);
    [U P Q] = lu_obj.U(pm_string);
  else
    % permutation vectors/matrices are not requested
    L = lu_obj.L0();
    U = lu_obj.U();
  end

end
