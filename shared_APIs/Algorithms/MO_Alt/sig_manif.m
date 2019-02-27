function [y, cost] = sig_manif(Fopt, FRF, FBB)
[Nt, NRF] = size(FRF);

global manifold;
problem.M = manifold;

f = Fopt(:);
A = kron(FBB.', eye(Nt));

problem.cost  = @(x) (f-A*x)'*(f-A*x);
problem.egrad = @(x) -2*A'*(f-A*x);

% checkgradient(problem);
warning('off', 'manopt:getHessian:approx');

[x,cost] = conjugategradient(problem,FRF(:));
% [x,cost,info,options] = trustregions(problem, FRF(:));
% info.iter
y = reshape(x,Nt,NRF);

end