function [ber,rate] = MMSE_method()

%the algorithm for fully_digital design
%cite the paper Generalized linear precoder and decoder design for MIMO channels using the
%weighted mmse criterion

 global Ns Vn H V_mopt W_mopt;
 [N_r,~] = size(H);
 Rnn = Vn *eye(N_r) ;
 S = H' * Rnn^(-1) * H ;
 [A,B] = eig(S);
 %% 特征值分解+降序排列取前N_RF列
 b = diag(B) ;    
 [B_sort, B_index] = sort(b,'descend');
 A = A(:,B_index);
 B = diag(B_sort);
 A = A(:,1:Ns);
 B = B(1:Ns,1:Ns);
 b = diag(B) ; 
 c = diag(B.^(-1/2));  %\lamda ^-1/2
 C = diag(c);
 d = diag(B.^(-1));
 D = diag(d);
 %% 求解 /mu + 1 
   n = Ns;
   while ( n>0) 
       u = sum(c(1:n))/(1+sum(d(1:n)));   %\mu^1/2
        if u <= b(n)
            break
        else 
            n = n - 1;
        end
   end
   
   f = (u^(-1) * C - D)^(0.5);
   M = zeros(Ns,Ns);
   M(1:n,1:n) = f(1:n,1:n);   %取前几列
   K = (u * C - u^2 * D )^(0.5) * C;
   Y = zeros(Ns,Ns);
   Y(1:n,1:n) = K(1:n,1:n);
   V_mopt = A * M;
   W_mopt = (Y * A' * H' * Rnn^(-1))';
   
   %% save results  
ber = get_ber(V_mopt,W_mopt);
rate = get_rate(V_mopt,W_mopt);