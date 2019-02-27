function BER = get_ber(V_equal, W_equal)

global Ns Nsym Nr Vn H hMod hDemod;
%one channel realization for Nsym times data streams transmission
% 2 is for real and imaginary
data = randi([0 1],Nsym*Ns*2,1);

%QPSK modulation generate all original signals s in Nsym times
s = reshape(step(hMod,data),Ns, Nsym);

%generate noise vector u
u = sqrt(Vn/2).*(randn(Nr,Nsym)+1i*randn(Nr,Nsym));

%get the receive vector 
%colloct Nsym receive vector in the r matrix, where nth column is a
%receive vector at n time
r = W_equal' * H * V_equal * s + W_equal' * u;

%QPSK demodulation
de_data = step(hDemod, r(:));
%get the number of error bits
ber = biterr(data,de_data);
BER = ber/length(data);
