function [L U P Q] = lusol(A,pm_opt,varargin)
  %LUSOL  performs LU factorization with LUSOL, returns L and U factors
  %
  % Use this function if you want direct access to the L and U factors as
  % computed by LUSOL fortran routines.
  %
  % Basic usage:
  %  [L U P Q] = lusol(A);
  %
  % By default P and Q are sparse permutation matrices.
  %
  % Get default options structure:
  %  lusol_options = lusol();
  %
  % Advanced usage:
  %  [L U p q] = lusol(A,'vector');
  %  [L U p q] = lusol(A,'vector',lusol_options);
  %  [L U p q] = lusol(A,'vector','pivot','TPP');
  %  [L U] = lusol(A);
  %
  % The second argument to lusol is a string specifying the desired form of the
  % permutation information.  Setting this to 'matrix' or [] tells lusol to
  % construct sparse permutation matrices.  This is the default behavior as it
  % most closely matches the behavior of other Matlab LU functions.  If
  % permutation vectors are desired, simply set this option to 'vector'.
  %
  % All subsequent options are LUSOL options.  These may be set in an options
  % structure or as key value pairs.
  %
  % If permutation variables are not specified, lusol will return permuted L and
  % U factors.
  %
  % See also: lusol_obj lusol_obj.luset
  %

  % returns options structure in special case
  if nargin == 0 && nargout <= 1
    L = lusol_obj.luset();
    return;
  end

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
