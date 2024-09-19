% 4 x 4 MIMO信道容量随信噪比的变化
% 采用极化天线与不采用极化天线的对比
clc;
clear;
U = 4;
S = 4;
I = eye(U,S);
snr = 1:20;
C = zeros(2,length(snr));
C_XPR = zeros(2,length(snr));
d2D = 250;
color = ['-o';'-s'];
color_XPR = [':o';':s'];
scene = 1;
    [sigma_ASD,sigma_ZSD,sigma_DS,sigma_SF,sigma_ASA,sigma_ZSA,N_cluster,N_ray,...
           c_ASD,m_ZSD,Dsp,Pcs,m_ZSD_offset,c_ZSA,c_ASA,m_xpr,s_xpr]=generate_para(scene,d2D);
    E_XPR = exp(s_xpr^2.*(log(10).^2)/200-m_xpr.*log(10)/10); 
    
    H = channel_3DMIMO_ULA(scene,d2D,U,S);
    H_MIMO = mean(H,3);
    
    H_XPR = channel_3DMIMO(scene,d2D,U,S);
    H_MIMO_XPR = mean(H_XPR,3);
    
    for SNR = 1:20 
        SNR_linear = 10.^(0.1*SNR/(1+E_XPR));
        SNR_linear_ULA = 10.^(0.1*SNR);
        C(scene,SNR) = log2(det(I+(SNR_linear_ULA/S)*(H_MIMO*H_MIMO.')));
        C(scene,SNR) = real(C(scene,SNR));
        C_XPR(scene,SNR) = log2(det(I+(SNR_linear/S)*(H_MIMO_XPR*H_MIMO_XPR.')));
        C_XPR(scene,SNR) = real(C_XPR(scene,SNR));
    end
plot(snr,C(1,:),'-o','LineWidth',2);
hold on;

plot(snr,C_XPR(1,:),'-^','LineWidth',2);
hold on;


C = zeros(2,length(snr));
C_XPR = zeros(2,length(snr));
scene = 2;
    [sigma_ASD,sigma_ZSD,sigma_DS,sigma_SF,sigma_ASA,sigma_ZSA,N_cluster,N_ray,...
           c_ASD,m_ZSD,Dsp,Pcs,m_ZSD_offset,c_ZSA,c_ASA,m_xpr,s_xpr]=generate_para(scene,d2D);
    E_XPR = exp(s_xpr^2.*(log(10).^2)/200-m_xpr.*log(10)/10); 
    
    H = channel_3DMIMO_ULA(scene,d2D,U,S);
    H_MIMO = mean(H,3);
    
    H_XPR = channel_3DMIMO(scene,d2D,U,S);
    H_MIMO_XPR = mean(H_XPR,3);
    
    for SNR = 1:20 
        SNR_linear = 10.^(0.1*SNR/(1+E_XPR));
        SNR_linear_ULA = 10.^(0.1*SNR);
        C(scene,SNR) = log2(det(I+(SNR_linear_ULA/S)*(H_MIMO*H_MIMO.')));
        C(scene,SNR) = real(C(scene,SNR));
        C_XPR(scene,SNR) = log2(det(I+(SNR_linear/S)*(H_MIMO_XPR*H_MIMO_XPR.')));
        C_XPR(scene,SNR) = real(C_XPR(scene,SNR));
    end
plot(snr,C_XPR(scene,:),':o','LineWidth',2);
hold on;
plot(snr,C(scene,:),':^','LineWidth',2);
hold on;



xlabel('SNR/(dB)');
ylabel('信道容量/(b/s/Hz)');
legend('UMi-ULA','UMi-XPR','UMa-ULA','UMa-XPR');
grid on;
