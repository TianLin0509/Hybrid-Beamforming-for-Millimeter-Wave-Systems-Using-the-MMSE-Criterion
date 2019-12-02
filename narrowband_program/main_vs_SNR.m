% Simulation codes for Hybrid Beamforming for Millimeter Wave Systems Using the MMSE Criterion
% Author : Lin Tian    2018. 08. 09 for paper revision
% this is corresponding to the fig: BER or rate v.s. SNR

clear all;  close all; clc;
disp(datestr(now));

%set up simulation parameters;
SNR_dB = (-15:5:5);

%All variables are corresponding to paper
%numbers of antennas, streams, RF chains, block
global Nt Nr Ns Nrf Nsym;
Nt = 64;
Nr = 64;
Ns =2;
Nrf = 2;
Nsym = 64;  %the number of symbols in one block


global N_loop;
N_loop = 500;   %iteration number

global Vn H Codebook_v Codebook_w n; %Vn £º Noise_power
%Codebook : the required codebook for OMP method

%using QPSK modulation
global hMod hDemod;
hMod = comm.PSKModulator(4,'BitInput',true,'PhaseOffset',pi/4);
hDemod = comm.PSKDemodulator('ModulationOrder',4,'BitOutput',true,'PhaseOffset',pi/4);
fprintf('params: \n Nt: %d  Nr: %d  Ns: %d N_loop: %d Nrf: %d \n SNR: %d : %d \n',...
    Nt,Nr,Ns,N_loop,Nrf,SNR_dB(1),SNR_dB(end));

%global initialization (optimal based on MMSE or rate)
global  V_mopt W_mopt V_ropt W_ropt ;
global manifold ;
%generate the manifold
manifold = complexcirclefactory(Nt*Nrf);

for snr_index = 1 : length(SNR_dB)
    Vn = 1 / 10^(SNR_dB(snr_index)/10);   % Noise Power    
    for  n = 1 : N_loop
        
        % generate channel matrix, codebooks for OMP
        [H ,Codebook_v, Codebook_w]  = channel_generation(Nt,Nr);
        %run different algorithms
        [Mrate_ber(snr_index,n), Mrate_rate(snr_index,n)] = Mrate_method();
        [MMSE_ber(snr_index,n), MMSE_rate(snr_index,n)] = MMSE_method();
        [Yuwei_ber(snr_index,n), Yuwei_rate(snr_index,n)] = Yuwei_method();
        [MO_Alt_ber(snr_index,n), MO_Alt_rate(snr_index,n)] = MO_Alt_method();
        [OMP_ber(snr_index,n), OMP_rate(snr_index,n)] = OMP_method();
        [MO_ber(snr_index,n), MO_rate(snr_index,n)] = MO_method();
        [GEVD_ber(snr_index,n), GEVD_rate(snr_index,n)] = GEVD_method();
        
    end   
    fprintf('current SNR: %d \n', SNR_dB(snr_index))
    disp(datestr(now));
end

%average Nloop channels
Mrate_Ber = mean(Mrate_ber, 2);  MMSE_Ber = mean(MMSE_ber, 2);  Yuwei_Ber = mean(Yuwei_ber, 2); MO_Alt_Ber = mean(MO_Alt_ber, 2);  OMP_Ber = mean(OMP_ber, 2); 
MO_Ber = mean(MO_ber, 2);   GEVD_Ber = mean(GEVD_ber, 2); 
Mrate_Rate = mean(Mrate_rate, 2);  MMSE_Rate = mean(MMSE_rate, 2); Yuwei_Rate = mean(Yuwei_rate, 2); MO_Alt_Rate = mean(MO_Alt_rate, 2); OMP_Rate = mean(OMP_rate, 2);
MO_Rate = mean(MO_rate, 2);  GEVD_Rate = mean(GEVD_rate, 2);
%plot figures for different metrics
figure(1)
plot(SNR_dB,  Mrate_Ber, 'm-pentagram ', 'LineWidth', 2)
hold on 
plot(SNR_dB,  MMSE_Ber, 'g-+ ', 'LineWidth', 2)
plot(SNR_dB,  Yuwei_Ber, 'b-o ', 'LineWidth', 2)
plot(SNR_dB,  MO_Alt_Ber, 'k-d ', 'LineWidth', 2)
plot(SNR_dB,  OMP_Ber, 'r-^ ', 'LineWidth', 2)
plot(SNR_dB,  MO_Ber, 'b--v ', 'LineWidth', 2)
plot(SNR_dB,  GEVD_Ber, 'k--x', 'LineWidth', 2)
%axis([SNR_dB(1), SNR_dB(end)])
xlabel('SNR(dB)')
ylabel('BER') 
legend('FD rate','FD MSE', 'HBF[19]',  'HBF[18]', 'OMP', 'MO', 'GEVD','northeast')
figure(2)
plot(SNR_dB, Mrate_Rate,'m-pentagram ', 'LineWidth', 2)
hold on 
plot(SNR_dB, MMSE_Rate,'g-+', 'LineWidth', 2)
plot(SNR_dB, Yuwei_Rate,'b-o', 'LineWidth', 2)
plot(SNR_dB, MO_Alt_Rate,'k-d', 'LineWidth', 2)
plot(SNR_dB,  OMP_Rate, 'r-^ ', 'LineWidth', 2)
plot(SNR_dB,  MO_Rate, 'b--v ', 'LineWidth', 2)
plot(SNR_dB,  GEVD_Rate, 'k--x', 'LineWidth', 2)
xlabel('SNR(dB)')
ylabel('Spectral Efficiency (bits/Hz/s)') 
legend('FD rate', 'FD MSE', 'HBF[19]', 'HBF[18]', 'OMP', 'MO', 'GEVD', 'Location', 'northwest')

