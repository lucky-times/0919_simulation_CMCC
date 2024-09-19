function NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id)
% 绘制UE上行发射功率的cdf
% 输入参数：
% SYS_config：LTE参数系统
% UEs：UE列表

for u_=1:length(activate_UE_id)
    UE_tx_power(:,u_) = UEs(activate_UE_id(u_)).ul_ad_TX_power; % 传入功控后的发射功率
end
UE_tx_power = reshape(UE_tx_power,1,[]);
UE_tx_power = 10*log10(UE_tx_power*1000); % 线性值转为dBm

if SYS_config.photo
    [range,output] = NR_plot_CDF(UE_tx_power);
    figure;
    plot(range,output);
    xlabel('UL Tx power(dBm)');
    ylabel('F(x)');
    title('UL ad frequency Tx power CDF');
end
%% 将功率数据保存为txt文档
UE_tx_power = sort(UE_tx_power);
index = 1:length(UE_tx_power)/100:length(UE_tx_power);
index = round(index);
R_ue_power = UE_tx_power(index);
R_ue_power = [R_ue_power max(UE_tx_power)];
fid=fopen('.\calibration\R_ue_power.txt','wt');
fprintf(fid,'%f\n',R_ue_power);
fprintf('该文档记录统计干扰UE功控后的发射功率\n');
fclose(fid);

end