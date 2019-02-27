function [ber, rate] = MO_method()

global  H  Vn W_mopt Nt  Nrf Nr;
i = 0;

W_equal = W_mopt;
w = trace (W_equal' * W_equal);
%random initialization
V_RF = exp( 1i*unifrnd(0,2*pi,Nt,Nrf));
W_RF = exp( 1i*unifrnd(0,2*pi,Nr,Nrf));
%iteration trigger, the normal initialization just for pass into functions
trigger = 1;
m_MSE_new = 100;

%limit the iterations number by i<10
while (trigger > 1e-5)
    
    % precoding
    H1 = H' * W_equal;
    Vn1 = Vn * w;
    [V_RF, V_U] = mo_algorithm(V_RF, Vn1, H1);
    V_equal = V_RF *V_U;
    v = trace (V_equal * V_equal');   %beta^(-2)
    
    %combining
    H2 = H * V_equal;
    Vn2 = Vn * v;
    [W_RF, W_B] = mo_algorithm(W_RF, Vn2, H2);
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



