function mse = get_mse(V_equal,W_equal,Weighted)
%compute mse or wmse

global Vn H Ns;

H_equal = W_equal' * H * V_equal;
E_matrix = H_equal * H_equal' - H_equal - H_equal' + eye(Ns) + Vn * W_equal'*W_equal;

%normal MSE
if nargin == 2
    mse = trace(E_matrix);
end

if nargin == 3
    mse = trace(Weighted * E_matrix);
end