function [V_RF, V_U, iter] = mo_algorithm(V_RF, Vn, H1)

global manifold;
[Nt, Nrf] = size(V_RF);
problem.M = manifold;

problem.cost = @(x)MMSE_cost(x,H1,Vn);
problem.egrad = @(x)MMSE_egrad(x,H1,Vn);

[x,iter] = conjugategradient(problem,V_RF(:));

V_RF = reshape(x,Nt,Nrf);
V_U = inv(V_RF'*H1 * H1'* V_RF+ 1 * Vn *(V_RF)'*V_RF)*V_RF'*H1;