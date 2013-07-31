classdef lusol_test < handle
  %LUSOL_TEST  test suite for lusol interface
  %
  % Calling lusol_test from the Matlab command line causes an lusol_test object
  % to be constructed.  The constructer calls all relevant test methods.
  %
  % Usage:
  %  >> lusol_test;
  %
  % Example output:
  %  >> lusol_test;
  %  test_addcol_lhr02: passed
  %  test_addrow_lhr02: passed
  %  test_delcol_lhr02: passed
  %  test_delrow_lhr02: passed
  %  test_factorize_lhr02: passed
  %  test_mulA_lhr02: passed
  %  test_mulAt_lhr02: passed
  %  test_r1mod_lhr02: passed
  %  test_repcol_lhr02: passed
  %  test_reprow_lhr02: passed
  %  test_solveA_lhr02: passed
  %  test_solveAt_lhr02: passed
  %  ... (more test methods may have been added)
  %

  properties

    data_dir = ''; % folder containing mat files

  end

  methods (Static)

    function tf = vectors_almost_equal(x,y,tol)
      %VECTORS_ALMOST_EQUAL  return true if input vectors are almost equal
      if nargin < 3 || isempty(tol)
        tol = 1e-6;
      end
      tf = norm(x-y,inf) <= tol;
    end

  end

  methods

    function obj = lusol_test()
      obj.set_data_dir();

      % get cell of method names
      method_cell = methods(obj);

      % loop through test methods
      for i = 1:length(method_cell)
        method_name = method_cell{i};
        if length(method_name) >= 5 && strcmp(method_name(1:5),'test_')
          fprintf([method_name ': ']);
          try
            tf = eval(['obj.' method_name]);
            if tf
              fprintf('passed\n');
            else
              fprintf('failed\n');
            end
          catch err
            fprintf('exception | %s | %s\n',err.identifier,err.message);
          end
        end
      end

    end

    function set_data_dir(obj)
      %SET_DATA_DIR  set appropriate data directory for test matrices
      s = which('lusol_test');
      obj.data_dir = strrep(s,'lusol_test.m','test_data/');
    end

    function A = load_matrix(obj,data_file)
      %LOAD_MATRIX  load a matrix from UF format .mat file
      S = load([obj.data_dir data_file]);
      A = S.Problem.A;
    end

    function tf = test_factorize_lhr02(obj)
      %TEST_FACTORIZE_LHR02  test LUSOL factorization of lhr02 matrix
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',100000);
      tf = mylu.inform() == 0;
    end

    function tf = test_solveA_lhr02(obj)
      %TEST_SOLVEA_LHR02  test LUSOL solve with lhr02 matrix
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',100000);
      [m n] = size(A);
      x = ones(n,1);
      b = A*x;
      xs = mylu.solveA(b);
      tf = obj.vectors_almost_equal(x,xs);
    end

    function tf = test_solveAt_lhr02(obj)
      %TEST_SOLVEA_LHR02  test LUSOL solve with lhr02 matrix transpose
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',100000);
      [m n] = size(A);
      x = ones(n,1);
      b = A'*x;
      xs = mylu.solveAt(b);
      tf = obj.vectors_almost_equal(x,xs);
    end

    function tf = test_mulA_lhr02(obj)
      %TEST_MULA_LHR02  test LUSOL multiply with lhr02 matrix
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',100000);
      [m n] = size(A);
      x = ones(n,1);
      b = A*x;
      bs = mylu.mulA(x);
      tf = obj.vectors_almost_equal(b,bs);
    end

    function tf = test_mulAt_lhr02(obj)
      %TEST_MULAT_LHR02  test LUSOL multiply with lhr02 matrix transpose
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',100000);
      [m n] = size(A);
      x = ones(n,1);
      b = A'*x;
      bs = mylu.mulAt(x);
      tf = obj.vectors_almost_equal(b,bs);
    end

    function tf = test_repcol_lhr02(obj)
      %TEST_REPCOL_LHR02  test LUSOL repcol with lhr02
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',100000);
      [m n] = size(A);
      x = ones(n,1);
      inform = mylu.repcol(x,1);
      tfv(1) = inform == 0;
      y = zeros(n,1);
      y(1) = 1;
      xs = mylu.mulA(y);
      tfv(2) = obj.vectors_almost_equal(x,xs);
      ys = mylu.solveA(x);
      tfv(3) = obj.vectors_almost_equal(y,ys);
      tf = all(tfv);
    end

    function tf = test_reprow_lhr02(obj)
      %TEST_REPROW_LHR02  test LUSOL reprow with lhr02
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',1000000);
      [m n] = size(A);
      x = ones(n,1);
      inform = mylu.reprow(x,1);
      tfv(1) = inform == 0;
      y = zeros(n,1);
      y(1) = 1;
      xs = mylu.mulAt(y);
      tfv(2) = obj.vectors_almost_equal(x,xs);
      ys = mylu.solveAt(x);
      tfv(3) = obj.vectors_almost_equal(y,ys);
      tf = all(tfv);
    end

    function tf = test_addcol_lhr02(obj)
      %TEST_ADDCOL_LHR02  test LUSOL addcol with lhr02
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',100000);
      [m nold] = size(A);
      % append a vector of ones
      x = ones(m,1);
      inform = mylu.addcol(x);
      tfv(1) = inform == 0;
      % check new size
      [m n] = mylu.size();
      tfv(2) = n == nold + 1;
      % check with a multiply
      y = zeros(n,1);
      y(n) = 1;
      xs = mylu.mulA(y);
      tfv(3) = obj.vectors_almost_equal(x,xs);
      tf = all(tfv);
    end

    function tf = test_addrow_lhr02(obj)
      %TEST_ADDROW_LHR02  test LUSOL addrow with lhr02
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',1000000);
      [mold n] = size(A);
      % append a vector of ones
      x = ones(1,n);
      inform = mylu.addrow(x);
      tfv(1) = inform == 0;
      % check new size
      [m n] = mylu.size();
      tfv(2) = m == mold + 1;
      % check with a multiply
      y = zeros(m,1);
      y(m) = 1;
      xs = mylu.mulAt(y);
      tfv(3) = obj.vectors_almost_equal(x(:),xs(:));
      tf = all(tfv);
    end

    function tf = test_delcol_lhr02(obj)
      %TEST_DELCOL_LHR02  test LUSOL delcol with lhr02
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',1000000);
      [m nold] = size(A);
      % delete column 1
      inform = mylu.delcol(1);
      tfv(1) = inform == -1;
      % check new size
      [m n] = mylu.size();
      tfv(2) = n == nold - 1;
      % check with a multiply
      y = ones(n,1);
      xs = mylu.mulA(y);
      x = A(:,2:end)*y;
      tfv(3) = obj.vectors_almost_equal(x(:),xs(:));
      tf = all(tfv);
    end

    function tf = test_delrow_lhr02(obj)
      %TEST_DELROW_LHR02  test LUSOL delrow with lhr02
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',1000000);
      [m n] = size(A);
      % delete row 1
      inform = mylu.delrow(1);
      tfv(1) = inform == -1;
      % check new size
      [m1 n1] = mylu.size();
      tfv(2) = (m1 == m) && (n1 == n);
      % check with a multiply
      y = zeros(m,1);
      y(end) = 1;
      xs = mylu.mulAt(y);
      x = zeros(n,1);
      tfv(3) = obj.vectors_almost_equal(x(:),xs(:));
      tf = all(tfv);
    end

    function tf = test_r1mod_lhr02(obj)
      %TEST_r1mod_LHR02  test LUSOL r1mod with lhr02
      A = obj.load_matrix('lhr02.mat');
      mylu = lusol_obj(A,'nzinit',1000000);
      [m n] = size(A);
      % set up rank 1 modification
      v = ones(m,1);
      w = ones(m,1);
      beta = 1.0;
      inform = mylu.r1mod(v,w,beta);
      tfv(1) = inform == 0;
      % compare with matlab matrix
      Ar = A + (beta*v)*w';
      b = ones(m,1);
      xm = Ar\b;
      % check with a solve
      x = mylu.solveA(b);
      tfv(2) = obj.vectors_almost_equal(x,xm);
      % check with multiply
      bs = mylu.mulA(xm);
      tfv(3) = obj.vectors_almost_equal(bs,b);
      tf = all(tfv);
    end

  end

end
