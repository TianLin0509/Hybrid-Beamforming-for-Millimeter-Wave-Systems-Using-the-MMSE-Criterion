function cost = MMSE_cost(x,H1,Vn)

global Nrf  Ns;
Nr = size(H1,1);

x = reshape(x,Nr,Nrf);

A = H1'*x ;
B = (x'*x);
C = B^(-1);
C = inv(B);
D = A*C;
E = D*x';
F = E*H1/Vn;
G = (F +eye(Ns))^(-1);

%cost = trace((H1'*x*(x'*x)^(-1)*x'*H1/Vn+eye(Ns))^(-1));
cost = trace(G);
end