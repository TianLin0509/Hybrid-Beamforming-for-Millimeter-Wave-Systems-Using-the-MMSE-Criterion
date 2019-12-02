function [H,AT,AR] = channel_generation(N_t,N_r)

%input: the numbers of transmit antennas and receive antennas
%output: the realized channel, codebook for vrf and wrf of OMP method 
%the comments is originally written in Chinese, hopefully I can have time
% to tranfer it to English, but not now.

N_c=5;
N_ray=10;
E_aoa = 2*pi* rand(N_c,1);                               %cluster的均值，服从（0,2*pi）的均匀分布
sigma_aoa = 10*pi/180;                                    %角度扩展为10°，化为弧度，即标准差
b = sigma_aoa/sqrt(2);                                    %根据标准差求相应的b，尺度参数
a = rand(N_c,N_ray)-0.5;                                  %生成(-0.5,0.5)区间内均匀分布的随机数列;
aoa = repmat(E_aoa,1,N_ray)-b*sign(a).*log(1-2*abs(a));   %生成符合拉普拉斯分布的随机数列(每列代表每个cluster的角度)
aoa = sin(aoa);
%-----------AOD
E_aod = 2*pi* rand(N_c, 1);                               %cluster的均值，服从（0,2*pi）的均匀分布
sigma_aod = 10*pi/180;                                    %角度扩展为10°，化为弧度，即标准差
b = sigma_aod/sqrt(2);                                    %根据标准差求相应的b，尺度参数
a = rand(N_c,N_ray)-0.5;                                  %生成(-0.5,0.5)区间内均匀分布的随机数列;
aod = repmat(E_aod,1, N_ray)-b*sign(a).*log(1-2*abs(a));   %生成符合拉普拉斯分布的随机数列(每行代表每个cluster的角度)
aod = sin(aod);


signature_t = [0:(N_t-1)]';
signature_t = 1i*pi* signature_t;                           %为接下来的signature构造做准备
signature_r = [0:(N_r-1)]';
signature_r = 1i*pi* signature_r;                           %为接下来的signature构造做准备

H_ray = zeros(N_r, N_t, N_c, N_ray);
H_cl = zeros(N_r, N_t, N_c);


    for i= 1: N_c
        for m = 1: N_ray
            H_ray(:,:,i,m)=complex(randn(1),randn(1))/sqrt(2)*exp((aoa(i,m)*signature_r))*exp((aod(i,m)*signature_t))'/sqrt(N_t*N_r); 
        end
    end  
        H_cl = sum(H_ray, 4);    
    
    H(:,:) = sqrt(N_t*N_r/N_c/N_ray)*sum(H_cl,3);   
    
    aod = aod(:).';
    aoa = aoa(:).';
    A = kron(aod,signature_t);
    AT = 1/sqrt(N_t)*exp(A);
    A = kron(aoa,signature_r);
    AR = 1/sqrt(N_r)*exp(A);
