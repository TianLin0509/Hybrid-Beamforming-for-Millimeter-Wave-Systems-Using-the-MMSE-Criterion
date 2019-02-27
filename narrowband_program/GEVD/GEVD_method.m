function [ber,rate] = GEVD_method()

%the proposed GEVD-HBF scheme

global H Vn  Nrf Nt Nr W_mopt;
i = 0;   %itertion index

%just to achieve W_RF*W_D = W_mopt
W_equal = W_mopt;
w = trace (W_equal' * W_equal);
%random initialization
V_RF = exp( 1i*unifrnd(0,2*pi,Nt,Nrf));
W_RF = exp( 1i*unifrnd(0,2*pi,Nr,Nrf));
%iteration trigger, the normal initialization just for pass into functions
trigger = 1;
m_MSE_new = 100;

%limit the iterations number by i<10
while (trigger > 1e-5 && i<10)
    
    % precoding
    H1 = H' * W_equal;
    [V_RF, V_U] = gevd_algorithm(V_RF, w, H1);
    V_equal = V_RF *V_U;
    v = trace (V_equal * V_equal');   %beta^(-2)
    
    %combining
    H2 = H * V_equal;
    [W_RF, W_B] = gevd_algorithm(W_RF, v, H2);
    W_equal = W_RF * W_B;
    w = trace (W_equal' * W_equal);
    %modified MSE
    H_equal = W_equal'*H2;
    
    m_MSE_old = m_MSE_new;
    m_MSE_new = trace(H_equal * H_equal' - H_equal - H_equal') + Vn * v * w;
    trigger = m_MSE_old - m_MSE_new;
    
    i = i + 1;
end

V_B = V_U / sqrt(v);

V = V_RF * V_B;
W = W_RF * W_B;

ber = get_ber(V, W);
rate = get_rate(V, W);



