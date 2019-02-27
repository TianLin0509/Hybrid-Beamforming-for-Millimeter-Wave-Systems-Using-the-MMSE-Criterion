function [V_RF,V_D] = MSEOMP (N_RF,H,AT,Vn1)
V_RF =[];
Ns = size(H,2);
VRES  = eye(Ns);
for i = 1 : N_RF
    vi = AT'*H*VRES;
    vi = diag(vi*vi');
    [a,k] = max(vi);
    V_RF =[V_RF AT(:,k)];
    AT(:,k) = [];
    V_D = inv(V_RF'*H*H'*V_RF+Vn1*(V_RF)'*V_RF)*V_RF'*H;
    RES = eye(Ns)-H'*V_RF*V_D;
    VRES = RES/norm(RES,'fro');
end
 
    
    
