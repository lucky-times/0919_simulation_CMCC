clear all;
clc
close all;
%% 加载UAV信息
str = '7.4.2多层窄波束配置--蜂窝三扇区组网/7.4.2多层窄波束配置--蜂窝三扇区组网-1km站间距-飞行高度500米-加扰50%.csv';
data = readtable(str);
info = data(:, [3,4,5,10,8]);
info = table2array(info);
temp = info((info(:, 1)>=480 & info(:, 1)<=525), :);% 高度、经度、纬度、SINR值

interCellPCI_measured = data(:, 11);
cellArray = table2array(interCellPCI_measured);
% 使用 cellfun 和 strsplit 来分割字符串
splittedStrings = cellfun(@(x) strsplit(x, ';'), cellArray, 'UniformOutput', false);
% 初始化一个空 cell 数组来存储数字
numbers = cell(size(cellArray));
% 遍历分割后的 cell 数组，将字符串转换为数字
for i = 1:numel(cellArray)
    for j = 1:numel(splittedStrings{i})
        numbers{i}(j) = str2double(splittedStrings{i}{j});
    end
end
interCellPCI_measured = numbers((info(:, 1)>=480 & info(:, 1)<=525), :);
cellSizes = cellfun(@numel, interCellPCI_measured);
% 找出最大的元素数目
maxSize = max(cellSizes);
cdfplot(cellSizes)



data2 = readtable('7个4.9G站点.xlsx');
info = data2(:, [9, 10]);
info = table2array(info);
pos = [];
for i = 1:length(info)
    if(ismember(info(i, 1), pos))
        continue;
    else
        pos = [pos;info(i, :)];
    end
end

%绘制UAV轨迹
figure(1)
scatter(temp(:,2), temp(:,3), 'r');
% 启用数据光标模式
dcm = datacursormode(gcf);
datacursormode on;
% 自定义数据光标显示的回调函数
set(dcm, 'UpdateFcn', @myupdatefcn);
% 显示图形
grid on; % 打开网格线
hold on;
% 绘制基站位置
scatter(pos(:,1), pos(:,2), 'b', 'filled');
xlabel('longitude')
ylabel('latitude')
legend('UAV','BS')
title('UAV h=500')
%计算UAV轨迹各个点到飞行圆心的距离
d = zeros(length(temp),1);
for i = 1:length(temp)
    d(i) = calculate_diatance(temp(i, 3), temp(i, 2));
end
hold off;

%% 绘制100m飞行高度，5种半径下的SINR值，100m,300m,500m,800m,1200m
h = 500;
r_100 = 31:91;
r_300 = 124:256;
r_500 = 288:492;
r_800 = 541:864;
r_1200 = 1090:1553;

% %% 处理PCI数据，查看是否有吻合的
% r100_interCellPCI_measured = interCellPCI_measured(r_100, :);
% writecell(r100_interCellPCI_measured, 'PCI_100m_Measured.csv');
% 
% r300_interCellPCI_measured = interCellPCI_measured(r_300, :);
% writecell(r300_interCellPCI_measured, 'PCI_300m_Measured.csv');
% 
% r500_interCellPCI_measured = interCellPCI_measured(r_500, :);
% writecell(r500_interCellPCI_measured, 'PCI_500m_Measured.csv');
% 
% r800_interCellPCI_measured = interCellPCI_measured(r_800, :);
% writecell(r800_interCellPCI_measured, 'PCI_800m_Measured.csv');
% 
% r1200_interCellPCI_measured = interCellPCI_measured(r_1200, :);
% writecell(r1200_interCellPCI_measured, 'PCI_1200m_Measured.csv');


%% 开始绘图
labelMode = 'Chinese'; %Chinese  English
linewidth1 = 2;
linewidth2 = 2;
axisLineWidth = 1;

r = 100;
figure
plot(1:length(r_100), temp(r_100, 4), 'b-', LineWidth=linewidth1);
hold on;
load("simulated_data_100.mat")
plot(1:length(plot_Sinr), plot_Sinr, 'r--', LineWidth=linewidth2);
xlim([1, length(plot_Sinr)+1]);
ylim([-17, 15]);
hold off;
diff = temp(r_100, 4)' - plot_Sinr;
mean1 = mean(diff);
set(gca, 'linewidth', axisLineWidth);
if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/(dB)')
    str = sprintf('Comparison of SINR for UAV at an Altitude of %dm and a Radius of %dm', h, r);
    legend('Measured', 'Simulated');
    str_filename = sprintf('Eng_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(gca, 'Position', [0.1, 0.12, 0.85, 0.8]); 
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('UAV在%dm高度，%dm半径下的SINR比较', h, r);
    legend('实测', '仿真'); 
    str_filename = sprintf('Chi_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 10, 'FontName', 'Microsoft YaHei UI');
    set(gca, 'Position', [0.1, 0.1, 0.85, 0.85]); 
end
title(str)
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3);
% 放大坐标区至充满图窗
axis tight;
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴
% 设置分辨率为300dpi
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', '-r300', str_filename);


r = 300;
figure
plot(1:length(r_300), temp(r_300, 4), 'b-', LineWidth=linewidth1);
hold on;
load("simulated_data_300.mat")
plot(1:length(plot_Sinr), plot_Sinr, 'r--', LineWidth=linewidth2);
xlim([1, length(plot_Sinr)+1]);
ylim([-10, 20]);
hold off;
diff = temp(r_300, 4)' - plot_Sinr;
mean2 = mean(diff);
set(gca, 'linewidth', axisLineWidth);
if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/(dB)')
    str = sprintf('Comparison of SINR for UAV at an Altitude of %dm and a Radius of %dm', h, r);
    legend('Measured', 'Simulated');
    str_filename = sprintf('Eng_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(gca, 'Position', [0.1, 0.12, 0.85, 0.8]); 
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('UAV在%dm高度，%dm半径下的SINR比较', h, r);
    legend('实测', '仿真'); 
    str_filename = sprintf('Chi_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 10, 'FontName', 'Microsoft YaHei UI');
    set(gca, 'Position', [0.1, 0.1, 0.85, 0.85]); 
end
title(str)
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3);
% 放大坐标区至充满图窗
axis tight;
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴
% 设置分辨率为300dpi
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', '-r300', str_filename);


r = 500;
figure
plot(1:length(r_500), temp(r_500, 4), 'b-', LineWidth=linewidth1);
hold on;
load("simulated_data_500.mat")
plot(1:length(plot_Sinr), plot_Sinr, 'r--', LineWidth=linewidth2);
xlim([1, length(plot_Sinr)+1]);
ylim([-10, 20]);
hold off;
diff = temp(r_500, 4)' - plot_Sinr;
mean3 = mean(diff);
set(gca, 'linewidth', axisLineWidth);
if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/(dB)')
    str = sprintf('Comparison of SINR for UAV at an Altitude of %dm and a Radius of %dm', h, r);
    legend('Measured', 'Simulated');
    str_filename = sprintf('Eng_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(gca, 'Position', [0.1, 0.12, 0.85, 0.8]); 
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('UAV在%dm高度，%dm半径下的SINR比较', h, r);
    legend('实测', '仿真'); 
    str_filename = sprintf('Chi_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 10, 'FontName', 'Microsoft YaHei UI');
    set(gca, 'Position', [0.1, 0.1, 0.85, 0.85]); 
end
title(str)
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3);
% 放大坐标区至充满图窗
axis tight;
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴
% 设置分辨率为300dpi
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', '-r300', str_filename);



r = 800;
figure
plot(1:length(r_800), temp(r_800, 4), 'b-', LineWidth=linewidth1);
hold on;
load("simulated_data_800.mat")
plot(1:length(plot_Sinr), plot_Sinr, 'r--', LineWidth=linewidth2);
xlim([1, length(plot_Sinr)+1]);
ylim([-15, 25]);
hold off;
diff = temp(r_800, 4)' - plot_Sinr;
mean4 = mean(diff);
set(gca, 'linewidth', axisLineWidth);
if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/(dB)')
    str = sprintf('Comparison of SINR for UAV at an Altitude of %dm and a Radius of %dm', h, r);
    legend('Measured', 'Simulated');
    str_filename = sprintf('Eng_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(gca, 'Position', [0.1, 0.12, 0.85, 0.8]); 
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('UAV在%dm高度，%dm半径下的SINR比较', h, r);
    legend('实测', '仿真'); 
    str_filename = sprintf('Chi_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 10, 'FontName', 'Microsoft YaHei UI');
    set(gca, 'Position', [0.1, 0.1, 0.85, 0.85]); 
end
title(str)
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3);
% 放大坐标区至充满图窗
axis tight;
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴
% 设置分辨率为300dpi
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', '-r300', str_filename);


r = 1200;
figure
plot(1:length(r_1200), temp(r_1200, 4), 'b-', LineWidth=linewidth1);
hold on;
load("simulated_data_1200.mat")
plot(1:length(plot_Sinr), plot_Sinr, 'r--', LineWidth=linewidth2);
xlim([1, length(plot_Sinr)+1]);
ylim([-15, 25]);
hold off;
diff = temp(r_1200, 4)' - plot_Sinr;
mean5 = mean(diff);
set(gca, 'linewidth', axisLineWidth);
if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/(dB)')
    str = sprintf('Comparison of SINR for UAV at an Altitude of %dm and a Radius of %dm', h, r);
    legend('Measured', 'Simulated');
    str_filename = sprintf('Eng_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman');
    set(gca, 'Position', [0.1, 0.12, 0.85, 0.8]); 
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('UAV在%dm高度，%dm半径下的SINR比较', h, r);
    legend('实测', '仿真'); 
    str_filename = sprintf('Chi_h%d_r%d.png', h, r);
    % 设置坐标轴字体
    set(gca, 'FontSize', 10, 'FontName', 'Microsoft YaHei UI');
    set(gca, 'Position', [0.1, 0.1, 0.85, 0.85]); 
end
title(str)
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3);
% 放大坐标区至充满图窗
axis tight;
% 设置分辨率为300dpi
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', '-r300', str_filename);
%%

figure
plot(1:length(r_100), temp(r_100, 4),LineWidth=2);
xlim([1,length(r_100)])
ylim([-8,15])

grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
set(gca, 'LineWidth', 1);
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴


%% 
if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/dB')
    str = sprintf('SINR Calculation Based on Field Data, h = 500m, r = 100m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 500m, r = 100m');
end  
title(str)

% figure
% plot(1:length(r_100), temp(r_100, 5),LineWidth=2.5);
% ylabel('PCI')
% title('r=100m')



figure
% subplot(5,1,2)
plot(1:length(r_300), temp(r_300, 4),LineWidth=2);
xlim([1,length(r_300)])
ylim([-10,20])

grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
set(gca, 'LineWidth', 1);
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴

if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/dB')
    str = sprintf('SINR Calculation Based on Field Data, h = 500m, r = 300m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 500m, r = 300m');
end  
title(str)

figure
% subplot(5,1,3)
plot(1:length(r_500), temp(r_500, 4),LineWidth=2);
xlim([1,length(r_500)])
ylim([-10,20]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
set(gca, 'LineWidth', 1);
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴

if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/dB')
    str = sprintf('SINR Calculation Based on Field Data, h = 500m, r = 500m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 500m, r = 500m');
end  
title(str)

figure
% subplot(5,1,4)
plot(1:length(r_800), temp(r_800, 4),LineWidth=2);
xlim([1,length(r_800)])
ylim([-15,25]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
set(gca, 'LineWidth', 1);
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴

if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/dB')
    str = sprintf('SINR Calculation Based on Field Data, h = 500m, r = 800m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 500m, r = 800m');
end  
title(str)

figure
% subplot(5,1,5)
plot(1:length(r_1200), temp(r_1200, 4),LineWidth=2);
xlim([1,length(r_1200)])
ylim([-15,25]);
grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
set(gca, 'LineWidth', 1);
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴

if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/dB')
    str = sprintf('SINR Calculation Based on Field Data, h = 500m, r = 1200m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 500m, r = 1200m');
end  
title(str)

save('real_Sinr.mat', "Sinr_real_r100","Sinr_real_r300","Sinr_real_r500","Sinr_real_r800","Sinr_real_r1200");

%% 回调函数定义
function txt = myupdatefcn(~, event_obj)
    % 获取点击点的索引
    pos = get(event_obj, 'Position');
    idx = get(event_obj, 'DataIndex');
    
    % 显示点的信息
    txt = {['X: ', num2str(pos(1))], ...
           ['Y: ', num2str(pos(2))], ...
           ['Index: ', num2str(idx)]};
end
