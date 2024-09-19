cd(SL_dir);
clear SNR I average_interval UEs
load('SISO_TU_SL.mat');
load('./data_files/SLvsLL_SISO_comparison.mat');
average_interval = UEs.walking_model.av_numb;
simulation_traces.UE_traces.TB_size = double(simulation_traces.UE_traces.TB_size);
for i = 1:average_interval:size(simulation_traces.UE_traces.SINR,3)-average_interval+1
    SNR(ceil(i/average_interval)) = simulation_traces.UE_traces.SNR(i);
    I1 = sum(double(simulation_traces.UE_traces.TB_size(1,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(1,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I2 = sum(double(simulation_traces.UE_traces.TB_size(2,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(2,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I(ceil(i/average_interval)) = I1+I2;
end
figure(1)
% clf
hold on
grid on
set(gca,'Fontsize',24);
plot(SNR,I,'ro-','Linewidth',2,'Markersize',14)
SNR_SL{1} = SNR;
I_SL{1} = I;

clear SNR I average_interval UEs
load('TxD_TU_SL.mat');
load('./data_files/SLvsLL_TxD_comparison.mat');
average_interval = UEs.walking_model.av_numb;
simulation_traces.UE_traces.TB_size = double(simulation_traces.UE_traces.TB_size);
for i = 1:average_interval:size(simulation_traces.UE_traces.SINR,3)-average_interval+1
    SNR(ceil(i/average_interval)) = simulation_traces.UE_traces.SNR(i);
    I1 = sum(double(simulation_traces.UE_traces.TB_size(1,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(1,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I2 = sum(double(simulation_traces.UE_traces.TB_size(2,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(2,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I(ceil(i/average_interval)) = I1+I2;
end
figure(1)
% clf
hold on
grid on
set(gca,'Fontsize',24);
plot(SNR,I,'bx-','Linewidth',2,'Markersize',16)
SNR_SL{2} = SNR;
I_SL{2} = I;

clear SNR I average_interval UEs
load('OLSM_TU_SL.mat');
load('./data_files/SLvsLL_OLSM_comparison.mat');
average_interval = UEs.walking_model.av_numb;
simulation_traces.UE_traces.TB_size = double(simulation_traces.UE_traces.TB_size);
for i = 1:average_interval:size(simulation_traces.UE_traces.SINR,3)-average_interval+1
    SNR(ceil(i/average_interval)) = simulation_traces.UE_traces.SNR(i);
    I1 = sum(double(simulation_traces.UE_traces.TB_size(1,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(1,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I2 = sum(double(simulation_traces.UE_traces.TB_size(2,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(2,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I(ceil(i/average_interval)) = I1+I2;
end
figure(1)
% clf
hold on
grid on
set(gca,'Fontsize',24);
plot(SNR,I,'+-','Linewidth',2,'Markersize',16,'Color',[0,0.75,0])
SNR_SL{3} = SNR;
I_SL{3} = I;

load('CLSM_TU_SL.mat');
load('./data_files/SLvsLL_CLSM_comparison.mat');
clear SNR I
average_interval = UEs.walking_model.av_numb;
simulation_traces.UE_traces.TB_size = double(simulation_traces.UE_traces.TB_size);
for i = 1:average_interval:size(simulation_traces.UE_traces.SINR,3)-average_interval+1
    SNR(ceil(i/average_interval)) = simulation_traces.UE_traces.SNR(i);
    I1 = sum(double(simulation_traces.UE_traces.TB_size(1,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(1,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I2 = sum(double(simulation_traces.UE_traces.TB_size(2,i:i+average_interval-1)).*simulation_traces.UE_traces.ACK(2,i:i+average_interval-1))/10^-3/10^6/average_interval;
    I(ceil(i/average_interval)) = I1+I2;
end
figure(1)
% clf
hold on
grid on
set(gca,'Fontsize',24);
plot(SNR,I,'d-','Linewidth',2,'Markersize',14,'Color',[0,0,0])
xlim([-10,40]);
SNR_SL{4} = SNR;
I_SL{4} = I;

cd(LL_dir)
load('SISO_TU_LL.mat');
plot(SNR_vec,mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6),'m--','Linewidth',2);
SNR_LL{1} = SNR_vec;
I_LL{1} = mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6);
load('TxD_TU_LL.mat');
plot(SNR_vec,mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6),'c--','Linewidth',2);
SNR_LL{2} = SNR_vec;
I_LL{2} = mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6);
load('OLSM_TU_LL.mat');
plot(SNR_vec,mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6),'g--','Linewidth',2);
SNR_LL{3} = SNR_vec;
I_LL{3} = mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6);
load('CLSM_TU_LL.mat');
plot(SNR_vec,mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6),'--','Linewidth',2,'Color',[0.5,0.5,0.5]);
SNR_LL{4} = SNR_vec;
I_LL{4} = mean((sum(simulation_results.cell_specific.throughput_coded,3))/LTE_params.Tsubframe/1e6);

xlabel('SNR [dB]');
ylabel('Throughput [Mbit/s]');
legend('SISO SL','TxD SL','OLSM SL','CLSM SL','SISO LL','TxD LL','OLSM LL','CLSM LL','Location','Northwest');
save('LLvsSL_data.mat','SNR_SL','SNR_LL','I_SL','I_LL');