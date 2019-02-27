function [ber,rate] = Mrate_method()
% traditional SVD algorithm for rate maximization

global  H Ns V_ropt W_ropt;
[U,~,V] = svd(H);
V_ropt = V(:,1:Ns);
%power constraint
V_ropt = V_ropt / norm(V_ropt,'fro');
W_ropt = U(:,1:Ns);

ber = get_ber(V_ropt,W_ropt);
rate = get_rate(V_ropt,W_ropt);