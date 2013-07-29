clear A
clear mylu
lusol_obj.unload()
A = rand(5);
mylu = lusol_obj(A);

%keyboard

% test various methods
[m n] = mylu.size();
stats = mylu.stats();
inform = mylu.inform();
nsing = mylu.nsing();
rank = mylu.rank();
depcol = mylu.depcol();
p = mylu.p();
q = mylu.q();
[U p q] = mylu.U();
[L0 p] = mylu.L0();

% test a simple solve
b = ones(5,1);
y = mylu.solveA(b);
A*y-b

% test a column replacement
[inform diag vnorm] = mylu.repcol(b,5);
Ar = [A(:,1:4) b];
yr = mylu.solveA(b);
Ar*yr - b