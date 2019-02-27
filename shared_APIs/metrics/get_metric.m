function obj = get_metric(obj)

%initialization
% the FD case
global Metric n;


 if obj.V_RF == 0
    V_equal = obj.V_B;
    W_equal = obj.W_B;

 else
    V_equal = obj.V_RF * obj.V_B;
    W_equal = obj.W_RF * obj.W_B;
end

    if (Metric.rate)
        obj.rate(n) = get_rate(V_equal, W_equal);
    end
    
    if (Metric.mse)
        obj.mse(n) = get_mse(V_equal, W_equal);
    end
    
     if (Metric.ber)
        obj.ber(n) = get_ber(V_equal, W_equal);
    end
    



