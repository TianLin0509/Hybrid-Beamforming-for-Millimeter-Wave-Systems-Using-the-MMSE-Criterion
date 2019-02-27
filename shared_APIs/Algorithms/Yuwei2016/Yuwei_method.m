function [ber, rate] = Yuwei_method()
% tHe YUWEI algoritHm for botH narrowband and broadband
% cite tHe paper Hybrid digital and analog beamforming design for large-scale antenna arrays
global Ns Vn Nrf Nt H ;

V_RF = yuweiA1();
Q = (V_RF'*V_RF);
T = Q^(-0.5);
L = H*V_RF*T;
[~,D,V] = svd(L);
[~,IX] = sort(diag(D),'descend');
M = V(:,IX);
U = M(:,1:Ns);
V_D = T*U;
V_D = V_D/norm(V_RF*V_D,'fro');

W_RF  = yuweiA2(V_D,V_RF);
J = W_RF'*H*V_RF*V_D*V_D'*V_RF'*H'*W_RF+Vn*W_RF'*W_RF;
W_D = J^(-1)*W_RF'*H*V_RF*V_D;


V = V_RF * V_D;
W = W_RF * W_D;

ber = get_ber(V, W);
rate = get_rate(V, W);
