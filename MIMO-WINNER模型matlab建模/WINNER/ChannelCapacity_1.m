% 不同规模MIMO信道容量对比
% 2 x 2 , 4 x 4 , 16 x 16 
clc;
clear;
U = [2,4,16];
S = [2,4,16];
snr = 1:20;
C_1 = zeros(3,length(snr));
C_2 = zeros(3,length(snr));
d2D = 250;

% UMi场景
scene = 1;
[sigma_ASD,sigma_ZSD,sigma_DS,sigma_SF,sigma_ASA,sigma_ZSA,N_cluster,N_ray,...
c_ASD,m_ZSD,Dsp,Pcs,m_ZSD_offset,c_ZSA,c_ASA,m_xpr,s_xpr]=generate_para(scene,d2D);
E_XPR = exp(s_xpr^2.*(log(10).^2)/200-m_xpr.*log(10)/10); 

for ii = 1:3
    I = eye(U(ii),S(ii));
    H = channel_3DMIMO(scene,d2D,U(ii),S(ii));
    H_MIMO = mean(H,3);
    
    for SNR = 1:20
        SNR_linear = 10.^(0.1*SNR/(1+E_XPR));
        C_1(ii,SNR) = log2(det(I+(SNR_linear/S(ii))*(H_MIMO*H_MIMO.')));
        C_1(ii,SNR) = real(C_1(ii,SNR));
    end
end

plot(snr,C_1(1,:),'-+','LineWidth',2);
hold on;
plot(snr,C_1(2,:),'-o','LineWidth',2);
hold on;
plot(snr,C_1(3,:),'-s','LineWidth',2);
hold on;
xlabel('SNR/(dB)');
ylabel('信道容量/(b/s/Hz)');
grid on;

% UMa场景
scene_2 = 2;
[sigma_ASD,sigma_ZSD,sigma_DS,sigma_SF,sigma_ASA,sigma_ZSA,N_cluster,N_ray,...
 c_ASD,m_ZSD,Dsp,Pcs,m_ZSD_offset,c_ZSA,c_ASA,m_xpr,s_xpr]=generate_para(scene_2,d2D);
E_XPR = exp(s_xpr^2.*(log(10).^2)/200-m_xpr.*log(10)/10); 

for ii = 1:3
    I = eye(U(ii),S(ii));
    H = channel_3DMIMO(scene_2,d2D,U(ii),S(ii));
    H_MIMO = mean(H,3);
    
    for SNR = 1:20
        SNR_linear = 10.^(0.1*SNR/(1+E_XPR));
        C_2(ii,SNR) = log2(det(I+(SNR_linear/S(ii))*(H_MIMO*H_MIMO.')));
        C_2(ii,SNR) = real(C_2(ii,SNR));
    end
end

plot(snr,C_2(1,:),':+','LineWidth',2);
hold on;
plot(snr,C_2(2,:),':o','LineWidth',2);
hold on;
plot(snr,C_2(3,:),':s','LineWidth',2);
legend('UMi-2x2','UMi-4x4','UMi-16x16','UMa-2x2','UMa-4x4','UMa-16x16');

