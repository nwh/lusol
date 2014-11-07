classdef lusol_obj < handle
  %LUSOL_OBJ  access to LUSOL for computations on a sparse LU factorization.
  %
  % Optional parameters:
  %  options = lusol_obj.luset();
  %  options = lusol_obj.luset('pivot','TRP');
  %
  % Initialization:
  %  mylu = lusol_obj(A);
  %
  % Initialization specifying options:
  %  mylu = lusol_obj(A,options);
  %
  % Note that initialization performs the factorization.  Subsequent
  % factorizations with the same object may be performed with the factorize
  % method.
  %
  % Factorize:
  %  [inform nsing depcol] = mylu.factorize(A);
  %  [inform nsing depcol] = mylu.factorize(A,options);
  %
  % There are several operations that may be performed with a lusol object.
  %
  % Usage:
  %  y = mylu.mulA(x);   % y = A*x
  %  y = mylu.mulAt(x);  % y = A'*x
  %  x = mylu.solveA(b); % x = A\b
  %  x = mylu.solveAt(b); % x = A'\b
  %
  % Update:
  %  inform = mylu.repcol(v,j); % replace column j of A with vector v
  %
  % See also:
  %  lusol
  %  lusol_obj.luset
  %  lusol_obj.factorize
  %  lusol_obj.stats
  %  lusol_obj.L0
  %  lusol_obj.U
  %  lusol_obj.p
  %  lusol_obj.q
  %  lusol_obj.mul
  %  lusol_obj.mulA
  %  lusol_obj.mulAt
  %  lusol_obj.mulL
  %  lusol_obj.mulLt
  %  lusol_obj.mulU
  %  lusol_obj.mulUt
  %  lusol_obj.solve
  %  lusol_obj.solveA
  %  lusol_obj.solveAt
  %  lusol_obj.solveL
  %  lusol_obj.solveLt
  %  lusol_obj.solveU
  %  lusol_obj.solveUt
  %  lusol_obj.repcol
  %  lusol_obj.reprow
  %  lusol_obj.addcol
  %  lusol_obj.addrow
  %  lusol_obj.delcol
  %  lusol_obj.delrow
  %  lusol_obj.r1mod
  %

  properties (Access=private)

    % object parameters

    nzinit = 0; % initial number of non zeros

    % lusol input parameters

    maxcol = 0; % max num cols searched for piv element (5)
    pivot = 0; % pivoting method [(0=TPP),1=TRP,2=TCP,3=TSP]
    keepLU = 0; % keep the nonzeros, if 0, permutations are computed (1)
    Ltol1 = 0; % max Lij allowed, default depends on luparm(6)
    Ltol2 = 0; % max Lij allowed during updates
    small = 0; % absolute tolerance for treating reals as zero (eps^0.8)
    Utol1 = 0; % absolute tol for flagging small diags of U (eps^0.67)
    Utol2 = 0; % rel tol for flagging small diags of U (eps^0.67)
    Uspace = 0; % (3.0)
    dens1 = 0; % (0.3)
    dens2 = 0; % (0.5)

    % lusol parameter vectors

    luparm_ptr = 0; % vector of integer parameters (input and output)
    parmlu_ptr = 0; % vector of double parameters (input and output)

    % scalars

    m_ptr = 0; % number of rows
    n_ptr = 0; % number of columns

    nelem_ptr = 0; % number of elements in original matrix
    nzmax_ptr = 0; % maximum storage allocated

    % vectors of lenth nzmax

    a_ptr = 0; % main storage array
    indc_ptr = 0; % row indices
    indr_ptr = 0; % column indices

    % vectors of length m

    p_ptr = 0; % row permutation
    lenr_ptr = 0;
    locr_ptr = 0;
    iqloc_ptr = 0;
    ipinv_ptr = 0;

    % vectors of length n

    q_ptr = 0; % column permutation
    lenc_ptr = 0;
    locc_ptr = 0;
    iploc_ptr = 0;
    iqinv_ptr = 0;

    % other

    depcol_lx = 0; % logical index indicating dependent columns
    int_class = 'int64'; % integer class used for integer arrays
    int_ptr_class = 'int64Ptr'; % integer class for libpointers

  end

  methods (Static)

    function options = luset(varargin)
      %LUSET  process input to create lusol options structure
      %
      % This method uses Matlab's inputParser class to handle lusol options.
      % It handles empty input, input structures, and key-value lists just like
      % many other Matlab programs.
      %
      % To obtain a structure with default settings use:
      %
      %  options = lusol.luset();
      %
      % To create an options structre that will use threshold rook pivoting,
      % use:
      %
      %  options = lusol.luset('pivot','TRP')
      %
      % The best reference for the lusol parameters is currently the comments
      % for the lu1fac subroutine in lusol1.f.
      %
      % |--------+----------+----------------------------------------------------|
      % | param  |  default | description                                        |
      % |--------+----------+----------------------------------------------------|
      %
      % lusol_obj options
      % |--------+----------+----------------------------------------------------|
      % | nzinit |        0 | minimum length for storage arrays                  |
      % |--------+----------+----------------------------------------------------|
      %
      % LUSOL integer parameters
      % |--------+----------+----------------------------------------------------|
      % | maxcol |        5 | max num cols searched for piv element              |
      % | pivot  |    'TPP' | pivoting method {'TPP','TRP','TCP','TSP'}          |
      % | keepLU |        1 | keep the nonzeros, if 0, permutations are computed |
      % |--------+----------+----------------------------------------------------|
      %
      % LUSOL real parameters
      % |--------+----------+----------------------------------------------------|
      % | Ltol1  |     10.0 | max Lij allowed, default depends on pivot method   |
      % | Ltol2  |     10.0 | max Lij allowed during updates                     |
      % | small  |  eps^0.8 | absolute tolerance for treating reals as zero      |
      % | Utol1  | eps^0.67 | absolute tol for flagging small diags of U         |
      % | Utol2  | eps^0.67 | rel tol for flagging small diags of U              |
      % | Uspace |      3.0 |                                                    |
      % | dens1  |      0.3 |                                                    |
      % | dens2  |      0.5 |                                                    |
      % |--------+----------+----------------------------------------------------|

      % get the input parser
      in_parse = inputParser;

      % storage parameters
      in_parse.addParamValue('nzinit',0,@(x) x>=0);

      % lusol integer parameters
      in_parse.addParamValue('maxcol',5,@(x) x>=0);
      in_parse.addParamValue('pivot','TPP',@(x) ismember(x,{'TPP','TRP','TCP','TSP'}));
      in_parse.addParamValue('keepLU',1,@(x) ismember(x,[0 1]));

      % lusol real parameters
      in_parse.addParamValue('Ltol1',10.0,@(x) x>=0.0);
      in_parse.addParamValue('Ltol2',10.0,@(x) x>=0.0);
      in_parse.addParamValue('small',eps^0.8,@(x) x>=0.0);
      in_parse.addParamValue('Utol1',eps^0.67,@(x) x>=0.0);
      in_parse.addParamValue('Utol2',eps^0.67,@(x) x>=0.0);
      in_parse.addParamValue('Uspace',3.0,@(x) x>=0.0);
      in_parse.addParamValue('dens1',0.3,@(x) x>=0.0);
      in_parse.addParamValue('dens2',0.5,@(x) x>=0.0);

      % parse the input
      in_parse.parse(varargin{:});

      % obtain the output
      options = in_parse.Results;

    end

    function r = get_range(p,i,j)
      %GET_RANGE  return a subset array from a libpointer array

      % get the datatype of the libpointer object
      pdatatype = p.DataType;
      % create a new pointer to p[i-1]
      pp = p + (i-1);
      % change the size of pp
      pp.setdatatype(pdatatype,j-i+1);
      % get the return array
      r = pp.Value;

    end

    function v = get_value(p,i)
      %GET_VALUE  return a single value from a libpointer array

      v = lusol_obj.get_range(p,i,i);

    end

    function load_library
      %LOAD_LIBRARY  Load the clusol shared library

      % do nothing if the library is already loaded
      if libisloaded('libclusol')
        return;
      end

      switch lower(computer)
        case 'glnxa64'
          % load clusol with linux prototype file
          loadlibrary('libclusol',@libclusol_proto_glnxa64);
        case 'maci64'
          % load clusol with mac prototype file
          loadlibrary('libclusol',@libclusol_proto_maci64);
        otherwise
          % the interface has not been implemented for other systems
          error('lusol_obj:load_library', ...
            'clusol library not implemented for this architecture')
      end

    end

    function unload_library
      %UNLOAD_LIBRARY  Unload the clusol shared library

      if libisloaded('libclusol')
        unloadlibrary('libclusol');
      end

    end

    function y = vector_process(x,xn)
      %VECTOR_PROCESS check and prepare vector for use in update routines
      %
      % Usage:
      %  y = lusol_obj.vector_process(x,xn);
      %
      % Input:
      %  x = input vector
      %  xn = expected length of x
      %
      % Output:
      %  y = output vector (double, full, and columnar)
      %
      % Error:
      %  The function will throw an error if x is not a vector of length xn
      %

      % check input vector
      if ~isvector(x) || length(x) ~= xn
        error('lusol_obj:vector_process','input is not a vector of correct length.');
      end

      % orient and densify
      y = double(full(x(:)));

    end

  end

  methods (Access=private)

    function parse_options(obj,varargin)
      %PARSE_OPTIONS  process the options structure to set parameters.

      % use luset to parse user input
      options = obj.luset(varargin{:});

      % storage parameters
      obj.nzinit = options.nzinit;

      % lusol integer parameters
      obj.maxcol = options.maxcol;
      obj.keepLU = options.keepLU;

      % lusol double parameters
      obj.Ltol1 = options.Ltol1;
      obj.Ltol2 = options.Ltol2;
      obj.small = options.small;
      obj.Utol1 = options.Utol1;
      obj.Utol2 = options.Utol2;
      obj.Uspace = options.Uspace;
      obj.dens1 = options.dens1;
      obj.dens2 = options.dens2;

      % set the pivoting strategy
      switch options.pivot
        case 'TPP'
          obj.pivot = 0;
        case 'TRP'
          obj.pivot = 1;
        case 'TCP'
          obj.pivot = 2;
        case 'TSP'
          obj.pivot = 3;
        otherwise
          error('lusol:set_options','Unkown pivot strategy.')
      end
    end

    function set_options(obj)
      %SET_OPTIONS  allocate and assign parameters to LUSOL arrays
      %
      % LUSOL stores input and output scalar parameters in two vectors:
      %
      %  luparm is an integer array of length 30
      %  parmlu is a double array of length 30
      %
      % this method sets the LUSOL input parameters in the correct location
      % for the fortran calls.
      %

      % allocate parameter vectors
      luparm = zeros(30,1,obj.int_class);
      parmlu = zeros(30,1,'double');

      % set parameter values
      luparm(2) = cast(-1,obj.int_class);
      luparm(3) = cast(obj.maxcol,obj.int_class);
      luparm(6) = cast(obj.pivot,obj.int_class);
      luparm(8) = cast(obj.keepLU,obj.int_class);
      parmlu(1) = double(obj.Ltol1);
      parmlu(2) = double(obj.Ltol2);
      parmlu(3) = double(obj.small);
      parmlu(4) = double(obj.Utol1);
      parmlu(5) = double(obj.Utol2);
      parmlu(6) = double(obj.Uspace);
      parmlu(7) = double(obj.dens1);
      parmlu(8) = double(obj.dens2);

      % allocate and initialize pointers
      obj.luparm_ptr = libpointer(obj.int_ptr_class,luparm);
      obj.parmlu_ptr = libpointer('doublePtr',parmlu);

    end

    function allocate_and_copy(obj,A)
      %ALLOCATE_AND_COPY  allocate LUSOL storage arrays
      %
      % LUSOL operates on many arrays.  This method allocates all of them
      % to an appropriate size.
      %

      % get information about A
      m = size(A,1);
      n = size(A,2);
      nelem = nnz(A);
      % set storage sizes
      nzmax = max([2*nelem 10*m 10*n 10000 obj.nzinit]);
      % vectors of length nzmax
      a = zeros(nzmax,1);
      indc = zeros(nzmax,1,obj.int_class);
      indr = zeros(nzmax,1,obj.int_class);
      % extract data from A for use in LUSOL
      [indc_tmp indr_tmp a_tmp] = find(A);
      indc(1:nelem) = cast(indc_tmp,obj.int_class);
      indr(1:nelem) = cast(indr_tmp,obj.int_class);
      a(1:nelem) = a_tmp;
      % vectors of length m
      p = zeros(m,1);
      lenr = zeros(m,1);
      locr = zeros(m,1);
      iqloc = zeros(m,1);
      ipinv = zeros(m,1);
      % vectors of length n
      q = zeros(n,1);
      lenc = zeros(n,1);
      locc = zeros(n,1);
      iploc = zeros(n,1);
      iqinv = zeros(n,1);

      %-- allocate and initialize libpointer "arrays" --%
      % integer scalars
      obj.m_ptr = libpointer(obj.int_ptr_class,m);
      obj.n_ptr = libpointer(obj.int_ptr_class,n);
      obj.nelem_ptr = libpointer(obj.int_ptr_class,nelem);
      obj.nzmax_ptr = libpointer(obj.int_ptr_class,nzmax);
      % vectors of length nzmax
      obj.a_ptr = libpointer('doublePtr',a);
      obj.indc_ptr = libpointer(obj.int_ptr_class,indc);
      obj.indr_ptr = libpointer(obj.int_ptr_class,indr);
      % vectors of length m
      obj.p_ptr = libpointer(obj.int_ptr_class,p);
      obj.lenr_ptr = libpointer(obj.int_ptr_class,lenr);
      obj.locr_ptr = libpointer(obj.int_ptr_class,locr);
      obj.iqloc_ptr = libpointer(obj.int_ptr_class,iqloc);
      obj.ipinv_ptr = libpointer(obj.int_ptr_class,ipinv);
      % vectors of length n
      obj.q_ptr = libpointer(obj.int_ptr_class,q);
      obj.lenc_ptr = libpointer(obj.int_ptr_class,lenc);
      obj.locc_ptr = libpointer(obj.int_ptr_class,locc);
      obj.iploc_ptr = libpointer(obj.int_ptr_class,iploc);
      obj.iqinv_ptr = libpointer(obj.int_ptr_class,iqinv);

    end

    function update_check(obj)
      %UPDATE_CHECK  throw an error if this method is called after updates
      %
      % Some methods should not be called after updates to a factorization.
      % This method checks if any updated have occuered and throws and
      % error if this is the case.
      %

      nupdat = lusol_obj.get_value(obj.luparm_ptr,15);
      if nupdat > 0
        error('lusol:post_update_call_error', ...
              'nsing and depcol cannot be called after updates.')
      end

    end

    function [x inform resid] = clu6sol(obj,b,mode)
      %CLU6SOL  call lu6sol to perform various solves with L and U factors.
      %
      % Performs data processing and direct call to clusol function for
      % solves.  Right hand side must be a vector.
      %
      % This is a private method and should only be used by class methods.
      % Users should call the various public interface methods.
      %
      % Usage:
      %  x = obj.clu6sol(b,mode)
      %
      % Input:
      %  b = right hand side vector
      %  mode = solution mode (see table below)
      %
      % Output:
      %  x = solution vector
      %  inform = status flag
      %  resid = 1-norm of residual
      %
      % Modes:
      %  1    x  solves   L x = b
      %  2    x  solves   L'x = b
      %  3    x  solves   U x = b
      %  4    x  solves   U'x = b
      %  5    x  solves   A x = b (default)
      %  6    x  solves   A'x = b
      %
      % inform flags:
      %  0 = successful solve
      %  1 = if U is singular, and residual is non-zero
      %

      % handle function options
      if nargin < 3
        mode = 5;
      end
      % make sure b is a vector
      if ~isvector(b)
        error('lusol:solve','b must be a vector.')
      end
      % orient b vector
      b = double(full(b(:)));
      lenb = length(b);
      % get size
      [m n] = obj.size();
      % allocate v and w vectors
      v = zeros(m,1,'double');
      w = zeros(n,1,'double');
      switch mode
        case 1
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 2
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 3
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 4
          if lenb ~= n, error('lusol:solve','b has incorrect size.'); end
          w = b;
        case 5
          if lenb ~= m, error('lusol:solve','b has incorrect size.'); end
          v = b;
        case 6
          if lenb ~= n, error('lusol:solve','b has incorrect size.'); end
          w = b;
        otherwise
          error('lusol:solve','unrecognized mode.')
      end

      % set up local libpointers for function call
      v_ptr = libpointer('doublePtr',v);
      w_ptr = libpointer('doublePtr',w);
      mode_ptr = libpointer(obj.int_ptr_class,mode);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);

      % call lusol routine
      calllib('libclusol','clu6sol', ...
        mode_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr);

      switch mode
        case 1
          x = v_ptr.Value;
        case 2
          x = v_ptr.Value;
        case 3
          x = w_ptr.Value;
        case 4
          x = v_ptr.Value;
        case 5
          x = w_ptr.Value;
        case 6
          x = v_ptr.Value;
      end

      inform = double(lusol_obj.get_value(obj.luparm_ptr,10));
      resid = double(lusol_obj.get_value(obj.parmlu_ptr,20));
    end

    function y = clu6mul(obj,x,mode)
      %CLU6MUL  call LUSOL to perform various multiplies with L and U factors.
      %
      % Performs data processing and direct call to clusol function to compute
      % a matrix vector multiply with the desired factor.
      %
      % This is a private method and should only be used by class methods.
      % Users should call the various public interface methods.
      %
      % Usage:
      %  y = obj.clu6mul(x,mode)
      %
      % Input:
      %  x = vector to multiply
      %  mode = multiply mode, see table below
      %
      % Output:
      %  y = product of desired factor and x
      %
      % Modes:
      %  1    y = L x
      %  2    y = L'x
      %  3    y = U x
      %  4    y = U'x
      %  5    y = A x (default)
      %  6    y = A'x
      %

      % handle optional input
      if nargin < 3
        mode = 5;
      end

      % check if x is a vector
      if ~isvector(x)
        error('lusol:mul','x must be a vector.');
      end

      % orient x vector
      x = double(full(x(:)));
      lenx = length(x);

      % get matrix size
      [m n] = obj.size();

      % set up temporary vectors
      v = zeros(m,1);
      w = zeros(n,1);

      switch mode
        case 1
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = x;
        case 2
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = x;
        case 3
          if lenx ~= n, error('lusol:mul','x has incorrect size.'); end
          w = x;
        case 4
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = w;
        case 5
          if lenx ~= n, error('lusol:mul','x has incorrect size.'); end
          w = x;
        case 6
          if lenx ~= m, error('lusol:mul','x has incorrect size.'); end
          v = x;
        otherwise
          error('lusol:mul','unrecognized mode.')
      end

      % set up local libpointers
      mode_ptr = libpointer(obj.int_ptr_class,mode);
      v_ptr = libpointer('doublePtr',v);
      w_ptr = libpointer('doublePtr',w);

      % call clu6mul
      calllib('libclusol','clu6mul', ...
        mode_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr);

      switch mode
        case 1
          y = v_ptr.Value;
        case 2
          y = v_ptr.Value;
        case 3
          y = v_ptr.Value;
        case 4
          y = w_ptr.Value;
        case 5
          y = v_ptr.Value;
        case 6
          y = w_ptr.Value;
      end

    end

  end

  methods

    % constructor and main factorize method

    function obj = lusol_obj(A,varargin)
      %LUSOL_OBJ  constructor for lusol object, factorize A
      %
      % Creates lusol object and factorizes A.
      %
      % Example:
      %   mylu = lusol(A);
      %
      % See also: LUSOL_OBJ.FACTORIZE
      %

      % load the shared library
      obj.load_library;

      % factorize the matrix
      obj.factorize(A,varargin{:});

    end

    function [inform nsing depcol] = factorize(obj,A,varargin)
      %FACTORIZE  perform lu factorization on A
      %
      % This method tells LUSOL to perform an LU factorization on A.
      %
      % Usage (after mylu object is initialized):
      %  [inform nsing depcol] = mylu.factorize(A)
      %  [inform nsing depcol] = mylu.factorize(A,options)
      %
      % Input:
      %  A = matrix to factorize
      %  options = options structure (optional)
      %
      % Output:
      %  inform = status flag
      %  nsing = esimate of the number of singularities
      %  depcol = logical index of dependent columns
      %
      % inform code:
      %  0 if the LU factors were obtained successfully.
      %  1 if U appears to be singular, as judged by lu6chk.
      %  3 if some index pair indc(l), indr(l) lies outside
      %    the matrix dimensions 1:m , 1:n.
      %  4 if some index pair indc(l), indr(l) duplicates
      %    another such pair.
      %  7 if the arrays a, indc, indr were not large enough.
      %    Their length "lena" should be increased to at least
      %    the value "minlen" given in luparm(13).
      %  8 if there was some other fatal error.  (Shouldn't happen!)
      %  9 if no diagonal pivot could be found with TSP or TDP.
      %    The matrix must not be sufficiently definite
      %    or quasi-definite.
      %

      % parse and set options
      obj.parse_options(varargin{:});
      obj.set_options();

      % allocate arrays and copy data from matrix A
      obj.allocate_and_copy(A);

      % get matrix size
      [m n] = obj.size();

      % temporary storage
      %cols_ptr = libpointer(obj.int_ptr_class,zeros(n,1));
      %markc_ptr = libpointer(obj.int_ptr_class,zeros(n,1));
      %markr_ptr = libpointer(obj.int_ptr_class,zeros(m,1));
      w_ptr = libpointer('doublePtr',zeros(n,1));

      % run lusol
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);
      calllib('libclusol','clu1fac', ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        obj.nelem_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        obj.iploc_ptr, ...
        obj.iqloc_ptr, ...
        obj.ipinv_ptr, ...
        obj.iqinv_ptr, ...
        w_ptr, ...
        ret_inform_ptr);

      % error checking
      ret_inform = ret_inform_ptr.Value;
      switch ret_inform
        case 0
          % ok, LU factors obtained
        case 1
          % ok, LU factors obtained, rank deficient
        case 3
          % error, some index pair indc(l), indr(l) lies outside
          % the matrix dimensions 1:m , 1:n.  Should not happen, because
          % matlab controls the input to lu1fac.
          err = MException('lusol:factorize', ...
                           'LUSOL reports improper input. inform = %d',ret_inform);
          throw(err);
        case 4
          % error, some index pair indc(l), indr(l) duplicates
          % another such pair.  Should not happen, because
          % matlab controls the input to lu1fac.
          err = MException('lusol:factorize', ...
                           'LUSOL reports improper input. inform = %d',ret_inform);
          throw(err);
        case 7
          % error, not enough storage.  User needs to increase nzinit
          % parameter
          err = MException('lusol:factorize', ...
                           ['LUSOL needs more storage. ' ...
                            'Increase the nzinit parameter.  inform = %d'], ...
                           ret_inform);
          throw(err);
        case 8
          % some other fatal error
          err = MException('lusol:factorize', ...
                           'LUSOL other fatal error. inform = %d',ret_inform);
          throw(err);
        case 9
          % error, no diagonal pivot could be found with TSP or TDP.
          % The matrix must not be sufficiently definite or quasi-definite
          err = MException('lusol:factorize', ...
                           ['LUSOL no diagonal pivot could be found with ' ...
                            'TSP or TDP. inform = %d'], ...
                           ret_inform);
          throw(err);
      end

      % compute logical index of dependent columns
      obj.depcol_lx = (w_ptr.Value <= 0.0);

      % user requests inform flag
      if nargout > 0
        inform = obj.inform();
      end

      % user requests number of singularities
      if nargout > 1
        nsing = obj.nsing();
      end

      % user requests dependent column indicator
      if nargout > 2
        depcol = obj.depcol_lx;
      end

    end

    % methods to collect information about matrix and factorization

    function [m n] = size(obj)
      %SIZE  get size of factorized matrix
      m = double(obj.m_ptr.Value);
      n = double(obj.n_ptr.Value);
    end

    function s = stats(obj)
      %STATS  return LUSOL stats structure
      %
      % This method builds a Matlab struct containing LUSOL output
      % parameters.
      %
      % LUSOL output parameters:
      %
      % inform   Return code from last call to any LU routine.
      % nsing    No. of singularities marked in the
      %          output array w(*).
      % jsing    Column index of last singularity.
      % minlen   Minimum recommended value for  lena.
      % maxlen   ?
      % nupdat   No. of updates performed by the lu8 routines.
      % nrank    No. of nonempty rows of U.
      % ndens1   No. of columns remaining when the density of
      %          the matrix being factorized reached dens1.
      % ndens2   No. of columns remaining when the density of
      %          the matrix being factorized reached dens2.
      % jumin    The column index associated with DUmin.
      % numL0    No. of columns in initial  L.
      % lenL0    Size of initial  L  (no. of nonzeros).
      % lenU0    Size of initial  U.
      % lenL     Size of current  L.
      % lenU     Size of current  U.
      % lrow     Length of row file.
      % ncp      No. of compressions of LU data structures.
      % mersum   lu1fac: sum of Markowitz merit counts.
      % nUtri    lu1fac: triangular rows in U.
      % nLtri    lu1fac: triangular rows in L.
      % Amax     Maximum element in  A.
      % Lmax     Maximum multiplier in current  L.
      % Umax     Maximum element in current  U.
      % DUmax    Maximum diagonal in  U.
      % DUmin    Minimum diagonal in  U.
      % Akmax    Maximum element generated at any stage
      %          during TCP factorization.
      % growth   TPP: Umax/Amax    TRP, TCP, TSP: Akmax/Amax.
      % resid    lu6sol: residual after solve with U or U'.
      %

      luparm = double(obj.luparm_ptr.Value);
      parmlu = double(obj.parmlu_ptr.Value);

      s.inform = luparm(10);
      s.nsing = luparm(11);
      s.jsing = luparm(12);
      s.minlen = luparm(13);
      s.maxlen = luparm(14);
      s.nupdat = luparm(15);
      s.nrank = luparm(16);
      s.ndens1 = luparm(17);
      s.ndens2 = luparm(18);
      s.jumin = luparm(19);
      s.numL0 = luparm(20);
      s.lenL0 = luparm(21);
      s.lenU0 = luparm(22);
      s.lenL = luparm(23);
      s.lenU = luparm(24);
      s.lrow = luparm(25);
      s.ncp = luparm(26);
      s.mersum = luparm(27);
      s.nUtri = luparm(28);
      s.nLtri = luparm(29);
      s.Amax = parmlu(10);
      s.Lmax = parmlu(11);
      s.Umax = parmlu(12);
      s.DUmax = parmlu(13);
      s.DUmin = parmlu(14);
      s.Akmax = parmlu(15);
      s.growth = parmlu(16);
      s.resid = parmlu(20);

    end

    function info = inform(obj)
      %INFORM  return code from last call to LUSOL routines

      %info = double(obj.luparm_ptr.Value(10));
      info = double(lusol_obj.get_value(obj.luparm_ptr,10));
    end

    function k = nsing(obj)
      %NSING  number of singularities marked in depcol
      %
      % This method may not be called after updates.  The nsing parameter
      % is only computed after a full factorize.
      %

      % this method only works if no updates have occured.
      obj.update_check();
      %k = double(obj.luparm_ptr.Value(11));
      k = double(lusol_obj.get_value(obj.luparm_ptr,11));
    end

    function k = rank(obj)
      %RANK  the rank of the matrix determined by the number of independent columns
      %
      % This method uses the LUSOL parameter nsing.  This is only computed
      % after a full factorize.  Thus this method should not be used to
      % determine the rank of a matrix after updates.  In that case look at
      % the flags that are returned by the update methods.
      %
      % Rank determination with LUSOL is more reliable under threshold rook
      % pivoting.  Example options:
      %
      %  options = lusol.luset('pivot','TRP','Ltol1',10)
      %
      n = double(obj.n_ptr.Value);
      nsing = obj.nsing();
      k = n-nsing;
    end

    function d = depcol(obj)
      %DEPCOL  logical vector indicating dependent columns
      %
      % This method may not be called after updates.  The method looks at
      % data that only relevant after a factorize.
      %

      % this method only works if no updated have occured.
      obj.update_check();
      d = obj.depcol_lx;
    end

    function ip = p(obj)
      %P  return row permutation vector
      m = double(obj.m_ptr.Value);
      ip = double(obj.p_ptr.Value(1:m));
    end

    function iq = q(obj)
      %Q  return column permutation vector
      n = double(obj.n_ptr.Value);
      iq = double(obj.q_ptr.Value(1:n));
    end

    % methods to get the matrix factors

    function [U p q] = U(obj,pm_opt)
      %U  get the upper triangular factor U
      %
      % Extract the U factor from LUSOL data and return as a Matlab sparse
      % matrix.
      %
      % Set up:
      %   mylu = lusol(A);
      %
      % Usage:
      %   U1 = mylu.U();
      %   [U2 p q] = mylu.U();
      %   [U2 P Q] = mylu.U('matrix');
      %
      % The first call returns U as a permuted triangle.  The second
      % returns U as an upper triangular matrix with permutation vectors p
      % and q.  The third call returns sparse permutation matrices P and Q.
      % The result of the three calls would produce:
      %   U1(p,q) == P*U1*Q == U2 == upper triangular
      %

      %
      % After a factorize (call to lu1fac) LUSOL stores U by rows at the
      % start of arrays a and indr.  lenr(1:m) stores the number of entries
      % in each row in original order.  locr(1:m) points to the beginning
      % of rows and is stored in original order.
      %
      % Special care must be taken when A is rank deficient.  LUSOL
      % actually stores lenU-nsing entries.  I suppose the extra nsing
      % contained in lenU could be for the zeros on the diagonal.  However,
      % LUSOL seems to handle these implicitly.
      %

      % permutation flag, set true if user desires upper triangular U and
      % permutation vectors
      permflg = false;
      % matrix flag, set true if user desires sparse permutation matrices
      % instead of vectors
      matrflg = false;
      % handle optional function input
      if nargout == 3
        permflg = true;
      end
      if nargin >= 2 && strcmpi(pm_opt,'matrix')
        matrflg = true;
      end
      % get basic matrix information
      [m n] = obj.size();
      p = obj.p();
      q = obj.q();
      s = obj.stats();
      % initialize arrays for U triplets
      ui = zeros(s.lenU-s.nsing,1,'double');
      uj = zeros(s.lenU-s.nsing,1,'double');
      ua = zeros(s.lenU-s.nsing,1,'double');
      % obtain required matrix data
      a = obj.a_ptr.Value;
      lenr = obj.lenr_ptr.Value;
      locr = obj.locr_ptr.Value;
      indr = obj.indr_ptr.Value;
      % array position pointers
      k1 = 1;
      k2 = 1;
      % loop through (rows?)
      for i = 1:s.nrank
        % get row index
        piv = p(i);
        % get length of row
        len = double(lenr(piv));
        % get location of row
        loc = double(locr(piv));
        % set end pointer
        k2 = k1+len-1;
        % load data into triplet arrays
        ui(k1:k2) = piv*ones(len,1);
        uj(k1:k2) = double(indr(loc:loc+len-1));
        ua(k1:k2) = double(a(loc:loc+len-1));
        % increment row start pointer
        k1 = k1+len;
      end
      % generate sparse matrix
      U = sparse(ui,uj,ua,m,n);
      % handle optional permutation of U
      if permflg
        % produce and return upper triangular U
        U = U(p,q);
      end
      % handle optional generation of permutation matrices
      if matrflg
        % construct and return sparse permutation matrices
        p = sparse((1:m)',p,1,m,m);
        q = sparse(q,(1:n)',1,n,n);
      end
    end

    function [L0 p] = L0(obj,pm_opt)
      %L0  get the initial lower triangular factor L0
      %
      % Extracts the initial lower triangular factor from the LUSOL data
      % structure.  LUSOL stores updates to L in product form, thus updates
      % are not included in L0.
      %
      % Set up:
      %   mylu.factorize(A);
      %
      % Usage:
      %   L1 = mylu.get_L0();
      %   [L2 p] = mylu.get_L0();
      %   [L2 P] = mylu.get_L0('matrix');
      %
      % The first call returns L1 as a permuted triangle.  The second and
      % third call will return L2 as a lower triangular matrix with
      % permutation vector p or matrix P.  The result will give:
      %   L1(p,p) == P*L1*P' == L2 == lower triangular
      %

      %
      % After a factorize (call to lu1fac) LUSOL stores non-trivial columns
      % of L at the end of a, indc, and indr.  lenc(1:numL0) stores the
      % number of entries in each column, not including the 1 on the
      % diagonal.  The negatives of the elements of L are stored in a.
      %
      % indc gives the row indices for non-zero elements
      % indr gives the column indices
      %
      % It turns out that off diagonal elements of L are stored in triplet
      % form at the end of a, indc, and indr.  It remains to add the ones
      % on the diagonal.
      %

      % permutation flag, set true if user desires upper triangular U and
      % permutation vectors
      permflg = false;
      % matrix flag, set true if user desires sparse permutation matrices
      % instead of vectors
      matrflg = false;
      % handle function options
      if nargout == 2
        permflg = true;
      end
      if nargin >= 2 && strcmpi(pm_opt,'matrix')
        matrflg = true;
      end
      % obtain information from object
      [m n] = obj.size();
      p = obj.p();
      s = obj.stats();
      lena = double(obj.nzmax_ptr.Value);
      % get matrix data
      a = obj.a_ptr.Value;
      indc = obj.indc_ptr.Value;
      indr = obj.indr_ptr.Value;
      % allocate arrays for triplet form
      li = zeros(s.lenL0+m,1);
      lj = zeros(s.lenL0+m,1);
      la = zeros(s.lenL0+m,1);
      % read the triplet form from LUSOL data
      li(1:s.lenL0) = double(indc(lena-s.lenL0+1:lena));
      lj(1:s.lenL0) = double(indr(lena-s.lenL0+1:lena));
      la(1:s.lenL0) = -double(a(lena-s.lenL0+1:lena));
      % add 1's along the diagonal
      li(s.lenL0+1:end) = (1:m)';
      lj(s.lenL0+1:end) = (1:m)';
      la(s.lenL0+1:end) = ones(m,1);
      % create matlab sparse matrix
      L0 = sparse(li,lj,la);
      % handle optional permutation
      if permflg
        % produce and return lower triangular L0
        L0 = L0(p,p);
      end
      % handle optional conversion of permuation vector to matrix
      if matrflg
        % construct and return sparse permutation matrix
        p = sparse((1:m)',p,1);
      end
    end

    % solve methods

    function [X inform resid] = solve(obj,B,mode)
      %SOLVE  solve systems with matrix factors
      %
      % This function solves all of the relavent systems of equations.  If
      % right hand side B is a matrix, it will solve for matrix X.
      %
      % Input:
      %  B = right hand size, can be a vector or a matrix
      %  mode = solution mode, see table below
      %
      % Output:
      %  X = solution matrix
      %  inform = status flag vector, one element for each column of B
      %  resid = 1-norm of residuals for each solve
      %
      % Modes:
      %  1    X  solves   L * X = B
      %  2    X  solves   L'* X = B
      %  3    X  solves   U * X = B
      %  4    X  solves   U'* X = B
      %  5    X  solves   A * X = B (default)
      %  6    X  solves   A'* X = B
      %
      % inform flags:
      %  0 = successful solve
      %  1 = if U is singular, and residual is non-zero
      %

      % set default mode
      if nargin < 3
        mode = 5;
      end

      % get size of factorized matrix
      [m n] = obj.size();

      % get size of B
      [Br Bc] = size(B);

      % compute size of X
      Xc = Bc;
      switch mode
        case 1 % X  solves   L X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        case 2 % X  solves   L'X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        case 3 % X  solves   U X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is n by Bc
          Xr = n;
        case 4 % X  solves   U'X = B
          % B must have n rows
          if Br ~= n, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        case 5 % X  solves   A X = B
          % B must have m rows
          if Br ~= m, error('lusol:solve','B has incorrect size.'); end
          % X is n by Bc
          Xr = n;
        case 6 % X  solves   A'X = B
          % B must have n rows
          if Br ~= n, error('lusol:solve','B has incorrect size.'); end
          % X is m by Bc
          Xr = m;
        otherwise
          error('lusol:solve','unrecognized mode.')
      end

      % allocate space for output X
      X = zeros(Xr,Xc);
      inform = zeros(1,Bc);
      resid = zeros(1,Bc);

      % compute solutions for all columns
      for j = 1:Bc
        [X(:,j) inform(j) resid(j)] = obj.clu6sol(B(:,j),mode);
      end

    end

    function [X inform resid] = solveA(obj,B)
      %SOLVEA  solve A*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,5);
    end

    function [X inform resid] = solveAt(obj,B)
      %SOLVEAT  solve A'*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,6);
    end

    function [X inform] = solveL(obj,B)
      %SOLVEL  solve L*X = B.
      %
      % See also: lusol.solve
      [X inform] = obj.solve(B,1);
    end

    function [X inform] = solveLt(obj,B)
      %SOLVELT  solve L'*X = B.
      %
      % See also: lusol.solve
      [X inform] = obj.solve(B,2);
    end

    function [X inform resid] = solveU(obj,B)
      %SOLVEU  solve U*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,3);
    end

    function [X inform resid] = solveUt(obj,B)
      %SOLVEUT  solve U'*X = B.
      %
      % See also: lusol.solve
      [X inform resid] = obj.solve(B,4);
    end

    % multiply methods

    function Y = mul(obj,X,mode)
      %MUL  compute matrix multiplies with various factors
      %
      % Perform matrix-vector or matrix-matrix multiply with desired matrix
      % factors.
      %
      % Usage:
      %  Y = mylu.mul(X,mode)
      %
      % Input:
      %  X = matrix or vector to compute multiply
      %  mode = multiply mode, see table below
      %
      % Output:
      %  Y = product of desired factor and X
      %
      % Mode:
      %  1    Y = L * X
      %  2    Y = L'* X
      %  3    Y = U * X
      %  4    Y = U'* X
      %  5    Y = A * X (default)
      %  6    Y = A'* X
      %

      % set default mode
      if nargin < 3
        mode = 5;
      end

      % get size of X
      [Xr Xc] = size(X);

      % get size of factored matrix
      [m n] = obj.size();

      % compute size of X
      Yc = Xc;
      switch mode
        case 1 %  Y = L X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 2 %  Y = L'X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 3 %  Y = U X
          % X must have n rows
          if Xr ~= n, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 4 %  Y = U'X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is n by Xc
          Yr = n;
        case 5 %  Y = A X
          % X must have n rows
          if Xr ~= n, error('lusol:mul','X has incorrect size.'); end
          % Y is m by Xc
          Yr = m;
        case 6 %  Y = A'X
          % X must have m rows
          if Xr ~= m, error('lusol:mul','X has incorrect size.'); end
          % Y is n by Xc
          Yr = n;
        otherwise
          error('lusol:mul','unrecognized mode.')
      end

      % allocate space for output Y
      Y = zeros(Yr,Yc);

      % compute solutions for all columns
      for j = 1:Xc
        [Y(:,j)] = obj.clu6mul(X(:,j),mode);
      end

    end

    function Y = mulA(obj,X)
      %MULA  compute Y = A*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,5);
    end

    function Y = mulAt(obj,X)
      %MULAT  compute Y = A'*X.
      %
      % Warning: this does not seem to work at the moment.
      %
      % See also: lusol.mul
      Y = obj.mul(X,6);
    end

    function Y = mulL(obj,X)
      %MULL  compute Y = L*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,1);
    end

    function Y = mulLt(obj,X)
      %MULLT  compute Y = L'*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,2);
    end

    function Y = mulU(obj,X)
      %MULU  compute Y = U*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,3);
    end

    function Y = mulUt(obj,X)
      %MULUT  compute Y = U'*X.
      %
      % See also: lusol.mul
      Y = obj.mul(X,4);
    end

    % update methods

    function [inform diag vnorm] = repcol(obj,v,j)
      %REPCOL  update LU factorization to replace a column
      %
      % Usage:
      %  [inform diag vnorm] = mylu.repcol(v,j)
      %
      % Inputs:
      %  v = new column
      %  j = column to replace
      %
      % Outputs:
      %  inform = status flag
      %  diag = ?
      %  vnorm = ?
      %
      % On exit:
      %  inform = -1  if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  2  if the update seemed to be unstable
      %               (diag much bigger than vnorm).
      %  inform =  7  if the update was not completed (lack of storage).
      %  inform =  8  if j is not between 1 and n.
      %

      % get matrix size
      [m n] = obj.size();

      % check and process input vector
      v = lusol_obj.vector_process(v,m);

      % check index of column to replace
      if j < 1 || j > n
        error('lusol:repcol','j must be between 1 and n.');
      end

      % create pointer to data
      v_ptr = libpointer('doublePtr',v);

      % temporary storage
      w_ptr = libpointer('doublePtr',zeros(n,1));

      % set up pointers to pass to function
      mode1_ptr = libpointer(obj.int_ptr_class,1);
      mode2_ptr = libpointer(obj.int_ptr_class,1);
      j_ptr = libpointer(obj.int_ptr_class,j);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);
      diag_ptr = libpointer('doublePtr',0);
      vnorm_ptr = libpointer('doublePtr',0);

      % call the library function
      calllib('libclusol','clu8rpc', ...
        mode1_ptr, ...
        mode2_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        j_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr, ...
        diag_ptr, ...
        vnorm_ptr);

      inform = obj.inform();
      diag = diag_ptr.Value;
      vnorm = vnorm_ptr.Value;
    end

    function inform = reprow(obj,w,i)
      %REPROW  update LU factorization to replace a row
      %
      % Usage:
      %  inform = mylu.reprow(w,i)
      %
      % Inputs:
      %  w = new row
      %  i = row to replace
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform = -1  if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %  inform =  8  if i is not between 1 and m.
      %

      % get matrix size
      [m n] = obj.size();

      % check and process input vector
      w = lusol_obj.vector_process(w,n);

      % prepare libpointers
      mode1_ptr = libpointer(obj.int_ptr_class,1);
      mode2_ptr = libpointer(obj.int_ptr_class,1);
      i_ptr = libpointer(obj.int_ptr_class,i);
      v_ptr = libpointer('doublePtr',zeros(m,1));
      w_ptr = libpointer('doublePtr',zeros(n,1));
      wnew_ptr = libpointer('doublePtr',w);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);

      % call lusol function
      calllib('libclusol','clu8rpr', ...
        mode1_ptr, ...
        mode2_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        i_ptr, ...
        v_ptr, ...
        w_ptr, ...
        wnew_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr);

      % get the return value
      inform = obj.get_value(obj.luparm_ptr,10);

    end

    function [inform diag vnorm] = addcol(obj,v)
      %ADDCOL  update LU factorization to add column to end of A
      %
      % Usage:
      %  [inform diag vnorm] = mylu.addcol(v)
      %
      % Input:
      %  v = vector of length m to added to the end of A
      %
      % Outputs:
      %  inform = status flag
      %  diag = ?
      %  vnorm = ?
      %
      % On exit:
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % Note that this will change n and the size of some workspace vectors.
      %

      % get current matrix size
      [m nold] = obj.size();
      n = nold + 1;

      % check and process input vector
      v = lusol_obj.vector_process(v,m);

      % copy and increase size of "row" arrays
      % these arrays have length nold and must be increased to length n

      % allocate vectors of length n
      q = zeros(n,1,obj.int_class);
      lenc = zeros(n,1,obj.int_class);
      locc = zeros(n,1,obj.int_class);
      iploc = zeros(n,1,obj.int_class);
      iqinv = zeros(n,1,obj.int_class);

      % copy data
      q(1:nold) = obj.q_ptr.Value;
      lenc(1:nold) = obj.lenc_ptr.Value;
      %locc(1:nold) = obj.locc_ptr.Value; % locc is assumed to be zero here?
      iploc(1:nold) = obj.iploc_ptr.Value;
      iqinv(1:nold) = obj.iqinv_ptr.Value;

      % set values at location n
      q(n) = n;
      lenc(n) = 0;
      locc(n) = 0;

      % create libpointers
      obj.q_ptr = libpointer(obj.int_ptr_class,q);
      obj.lenc_ptr = libpointer(obj.int_ptr_class,lenc);
      obj.locc_ptr = libpointer(obj.int_ptr_class,locc);
      obj.iploc_ptr = libpointer(obj.int_ptr_class,iploc);
      obj.iqinv_ptr = libpointer(obj.int_ptr_class,iqinv);

      % set n_ptr
      obj.n_ptr = libpointer(obj.int_ptr_class,n);

      % prepare temporary data
      mode_ptr = libpointer(obj.int_ptr_class,1);
      v_ptr = libpointer('doublePtr',v);
      w_ptr = libpointer('doublePtr',zeros(n,1));
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);
      diag_ptr = libpointer('doublePtr',0);
      vnorm_ptr = libpointer('doublePtr',0);

      % call the clusol function
      calllib('libclusol','clu8adc', ...
        mode_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr, ...
        diag_ptr, ...
        vnorm_ptr);

      % prepare output
      inform = double(ret_inform_ptr.Value);
      diag = diag_ptr.Value;
      vnorm = vnorm_ptr.Value;

    end

    function [inform diag] = addrow(obj,w)
      %ADDROW  update LU factorization to add row to end of A
      %
      % Usage:
      %  [inform diag] = mylu.addrow(w)
      %
      % Input:
      %  w = vector of length n to added to the bottom of A
      %
      % Outputs:
      %  inform = status flag
      %  diag = ?
      %
      % On exit:
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % This will change m and the size of some workspace vectors.
      %

      % get matrix size
      [mold n] = obj.size();
      m = mold + 1;

      % check and process input vector
      w = lusol_obj.vector_process(w,n);

      % allocate vectors of length m
      p = zeros(m,1,obj.int_class);
      lenr = zeros(m,1,obj.int_class);
      locr = zeros(m,1,obj.int_class);
      iqloc = zeros(m,1,obj.int_class);
      ipinv = zeros(m,1,obj.int_class);

      % copy data
      p(1:mold) = obj.p_ptr.Value;
      lenr(1:mold) = obj.lenr_ptr.Value;
      locr(1:mold) = obj.locr_ptr.Value;
      iqloc(1:mold) = obj.iqloc_ptr.Value;
      ipinv(1:mold) = obj.ipinv_ptr.Value;

      % create libpointers
      obj.p_ptr = libpointer(obj.int_ptr_class,p);
      obj.lenr_ptr = libpointer(obj.int_ptr_class,lenr);
      obj.locr_ptr = libpointer(obj.int_ptr_class,locr);
      obj.iqloc_ptr = libpointer(obj.int_ptr_class,iqloc);
      obj.ipinv_ptr = libpointer(obj.int_ptr_class,ipinv);

      % set m_ptr
      obj.m_ptr = libpointer(obj.int_ptr_class,m);

      % prepare temporary data
      w_ptr = libpointer('doublePtr',w);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);
      diag_ptr = libpointer('doublePtr',0);

      % call the library function
      calllib('libclusol','clu8adr', ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr, ...
        diag_ptr);

      % gather the output
      inform = ret_inform_ptr.Value;
      diag = diag_ptr.Value;

    end

    function inform = delcol(obj,j)
      %DELCOL  update LU factorization to delete column j from A
      %
      % Usage:
      %  inform = mylu.delcol(j)
      %
      % Input:
      %  j = column to delete from A
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform =  -1 if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % This will decrement n by 1 and change relevant storage vectors.
      %

      % get matrix size
      [m n] = obj.size();

      % check to make sure there is more than 1 column
      if n == 1
        error('lusol:delcol','matrix has only one column, cannot delete.')
      end

      % check if j is within bounds
      if j < 1 || j > n
        error('lusol:delcol','j is out of bounds.');
      end

      % set up local data
      j_ptr = libpointer(obj.int_ptr_class,j);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);

      % call the lusol function
      calllib('libclusol','clu8dlc', ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        j_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr);

      % decrement n
      obj.n_ptr = libpointer(obj.int_ptr_class,n-1);

      % gather output
      inform = ret_inform_ptr.Value;

    end

    function inform = delrow(obj,i)
      %DELROW  update LU factorization to delete row i from A
      %
      % DELROW updates the LU factorization  A = L*U  when row  idel
      % (the vector  w) is deleted from  A. The update is implemented as
      % the rank-one modification
      %
      %           A(new)  =  A  -  e(idel) * w',
      %
      % followed by a renumbering that makes row  idel  the last row of  A
      % and shifts rows  idel + 1,  idel + 2,  ...,  m  one place up.
      % Thus, row  idel  is replaced by the zero vector (rather than being
      % deleted), and is then cyclically permuted to the bottom of  A.
      % The dimensions of  A  do not alter, but  A  and  U  may become
      % singular.
      %
      % Note --- significant overhead is involved in permuting row  idel
      % to the bottom.  In some cases it may be better to use  lu8rpr  to
      % replace row  idel  by zero, leaving it in the current position.
      % The growth of nonzeros in  L  and  U  is identical, but less
      % housekeeping is required than with  lu8dlr.
      %
      % Usage:
      %  inform = mylu.delrow(i)
      %
      % Input:
      %  i = row to delete from A
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform =  -1 if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %
      % This does not change the size of A and puts zeros in the last row.
      %

      % get the size of the matrix
      [m n] = obj.size();

      % check to make sure there is more than 1 row
      if m == 1
        error('lusol:delrow','matrix has only one row, cannot delete.')
      end

      % check input
      if i < 1 || i > m
        error(lusol:delrow,'i is out of bounds.');
      end

      % prepare local pointer data
      mode_ptr = libpointer(obj.int_ptr_class,1);
      i_ptr = libpointer(obj.int_ptr_class,i);
      v_ptr = libpointer('doublePtr',zeros(m,1));
      w_ptr = libpointer('doublePtr',zeros(n,1));
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);

      % call the library function
      calllib('libclusol','clu8dlr', ...
        mode_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        i_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr);

      % from LUSOL code comments.  This method does not change the size of
      % A.  It introduces a zero row, then permutes it to the bottom.  This
      % will have to be tested.  For now, I am not going to change the
      % value of m.

      inform = ret_inform_ptr.Value;

    end

    function inform = r1mod(obj,v,w,beta)
      %R1MOD  update LU factorization to perform rank-1 update A+beta*v*w'
      %
      % Usage:
      %  inform = mylu.r1mod(v,w,beta)
      %
      % Input:
      %  v = vector of length m
      %  w = vector of length n
      %  beta = scalar value
      %
      % Outputs:
      %  inform = status flag
      %
      % On exit:
      %  inform =  -1 if the rank of U decreased by 1.
      %  inform =  0  if the rank of U stayed the same.
      %  inform =  1  if the rank of U increased by 1.
      %  inform =  7  if the update was not completed (lack of storage).
      %

      % get size of matrix
      [m n] = obj.size();

      % handle variable input
      if nargin < 4 || isempty(beta)
        beta = 1.0;
      end

      % check and process input vectors
      v = lusol_obj.vector_process(v,m);
      w = lusol_obj.vector_process(w,n);

      % prepare temporary libpointers
      mode_ptr = libpointer(obj.int_ptr_class,1);
      beta_ptr = libpointer('doublePtr',beta);
      v_ptr = libpointer('doublePtr',v);
      w_ptr = libpointer('doublePtr',w);
      ret_inform_ptr = libpointer(obj.int_ptr_class,0);

      % call library function
      calllib('libclusol','clu8mod', ...
        mode_ptr, ...
        obj.m_ptr, ...
        obj.n_ptr, ...
        beta_ptr, ...
        v_ptr, ...
        w_ptr, ...
        obj.nzmax_ptr, ...
        obj.luparm_ptr, ...
        obj.parmlu_ptr, ...
        obj.a_ptr, ...
        obj.indc_ptr, ...
        obj.indr_ptr, ...
        obj.p_ptr, ...
        obj.q_ptr, ...
        obj.lenc_ptr, ...
        obj.lenr_ptr, ...
        obj.locc_ptr, ...
        obj.locr_ptr, ...
        ret_inform_ptr);

      % prepare function output
      inform = ret_inform_ptr.Value;

    end

    % inherited, not implemented methods

    function TF = eq(obj,x)
      %EQ  Not implemented, returns false
      TF = false;
    end

    function TF = ge(obj,x)
      %GE  Not implemented, returns false
      TF = false;
    end

    function TF = gt(obj,x)
      %GT  Not implemented, returns false
      TF = false;
    end

    function TF = le(obj,x)
      %LE  Not implemented, returns false
      TF = false;
    end

    function TF = lt(obj,x)
      %LT  Not implemented, returns false
      TF = false;
    end

    function TF = ne(obj,x)
      %NE  Not implemented, returns false
      TF = false;
    end

    function addlistener(obj,varargin)
      %ADDLISTENER  Not implemented
    end

    function notify(obj,varargin)
      %NOTIFY  Not implemented
    end

  end

end % classdef
