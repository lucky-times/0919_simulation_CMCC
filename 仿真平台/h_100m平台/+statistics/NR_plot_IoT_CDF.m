function NR_plot_IoT_CDF(SYS_config,IoT)
% 绘制IoTcdf曲线
% 输入参数：
% SYS_config：参数系统
% IoT：干扰对噪声的比值

if SYS_config.photo
    % 绘制IoT cdf曲线
    % 单系统的ToT中的干扰仅为同频干扰
    % 双系统中的干扰包含同频及邻频干扰
    IoT_1 = IoT;
    [range,output] = NR_plot_CDF(IoT_1);
    figure;
    plot(range,output,'b');
    grid on;
    xlabel('IoT（dB）');
    ylabel('F(x)');
    title('IoT CDF');
    % 保存结果，为了将不同功控敏感度的结果画在同一个fig上
    file_name1 = ['.\tmp\','IOT_range_',num2str((SYS_config.PCe_param)^2)];% 保存在tmp中，如果需要画不同功控错误情况下的图形，可以写个脚本调用
    file_name2 = ['.\tmp\','IOT_output_',num2str((SYS_config.PCe_param)^2)];
    save(file_name1,'range');
    save(file_name2,'output');
end

end