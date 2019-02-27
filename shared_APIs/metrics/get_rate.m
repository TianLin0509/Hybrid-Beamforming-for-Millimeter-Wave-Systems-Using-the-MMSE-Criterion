function rate = get_rate(V_equal, W_equal)
%get the rate (SE) for equivalent V and W
global Vn H Ns;
rate = log2(det(eye(Ns) + 1/Vn * pinv(W_equal) * H * V_equal * V_equal' * H' *W_equal));
%log2(det(eye(Ns) + 1/Vn * pinv(WC_opt(:,:,i)) * H(:,:,i) * VC_opt(:,:,i) * VC_opt(:,:,i)' * H(:,:,i)' * WC_opt(:,:,i)));
