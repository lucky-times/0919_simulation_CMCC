function NR_plot_SINR_CDF_half(UE_SINR_sum_dl,UE_SINR_sum_ul)%上行数据暂时不输出
[range,output] = NR_plot_CDF(UE_SINR_sum_dl);
% if SYS_config.photo
if true
    figure;
    plot(range,output,'b');
    grid on;
    xlabel('SINR（dB）');
    ylabel('F(x)');
    title('50% DL+50% UL_DL SINR CDF');
    grid on
end
%% 保存下行SINR数据
UE_SINR = sort(UE_SINR_sum_dl);
index = 1:length(UE_SINR)/100:length(UE_SINR);
index = round(index);
R_sinr = UE_SINR(index);
R_sinr = [R_sinr max(UE_SINR)];
fid=fopen('.\calibration\R_sinr_half.txt','wt');
fprintf(fid,'%f\n',R_sinr);
fclose(fid);

[range,output] = NR_plot_CDF(UE_SINR_sum_ul);
% if SYS_config.photo
if true
    subplot(1,2,2);
    plot(range,output);
    xlabel('SINR（dB）');
    ylabel('F(x)');
    title('50% UL+50% DL_UL SINR CDF');
    grid on
end

ul_UE_SINR = sort(UE_SINR_sum_ul);
index = 1:length(ul_UE_SINR)/100:length(ul_UE_SINR);
index = round(index);
R_sinr = ul_UE_SINR(index);
R_sinr = [R_sinr max(ul_UE_SINR)];
fid=fopen('.\calibration\R_sinr_ul_half.txt','wt');
fprintf(fid,'%f\n',R_sinr);
fclose(fid);
end