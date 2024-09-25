clear all;
clc
close all;
%% 加载UAV信息
str = '7.4.2多层窄波束配置--蜂窝三扇区组网/7.4.2多层窄波束配置--蜂窝三扇区组网-1km站间距-飞行高度200米-加扰50%.csv';
data = readtable(str);
info = data(:, [3,4,5,10,8]);
info = table2array(info);
temp = info((info(:, 1)>=190 & info(:, 1)<=215), :);% 高度、经度、纬度、SINR值

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
interCellPCI_measured = numbers((info(:, 1)>=190 & info(:, 1)<=215), :);
cellSizes = cellfun(@numel, interCellPCI_measured);
% 找出最大的元素数目
maxSize = max(cellSizes);
cdfplot(cellSizes)

data2 = readtable('7个4.9G站点.xlsx');
info2 = data2(:, [9, 10]);
info2 = table2array(info2);
pos = [];
for i = 1:length(info2)
    if(ismember(info2(i, 1), pos))
        continue;
    else
        pos = [pos;info2(i, :)];
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
title('UAV h=200')
%计算UAV轨迹各个点到飞行圆心的距离
d = zeros(length(temp),1);
for i = 1:length(temp)
    d(i) = calculate_diatance(temp(i, 3), temp(i, 2));
end
hold off;

%% 绘制h飞行高度，5种半径下的SINR值，100m,300m,500m,800m,1200m
h = 200;
r_100 = 92:148;
r_300 = 206:339;
r_500 = 377:604;
r_800 = 680:1045;
r_1200 = 1251:1804;



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
a1 = mean(temp(r_100, 4));
b1 = mean(plot_Sinr);
xlim([1, length(plot_Sinr)+1]);
ylim([-17, 15]);
hold off;
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
a1 = mean(temp(r_300, 4));
b1 = mean(plot_Sinr);
xlim([1, length(plot_Sinr)+1]);
ylim([-10, 20]);
hold off;
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
a1 = mean(temp(r_500, 4));
b1 = mean(plot_Sinr);
xlim([1, length(plot_Sinr)+1]);
ylim([-10, 20]);
hold off;
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
a1 = mean(temp(r_800, 4));
b1 = mean(plot_Sinr);
xlim([1, length(plot_Sinr)+1]);
ylim([-15, 25]);
hold off;
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
a5 = mean(temp(r_1200, 4));
b5 = mean(plot_Sinr);
xlim([1, length(plot_Sinr)+1]);
ylim([-15, 25]);
hold off;
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

%% 求出加权因子
Sinr_simul_r100 = load('Sinr_simul_r100.mat');
Sinr_simul_r300 = load('Sinr_simul_r300.mat');
Sinr_simul_r500 = load('Sinr_simul_r500.mat');
Sinr_simul_r800 = load('Sinr_simul_r800.mat');
Sinr_simul_r1200 = load('Sinr_simul_r1200.mat');
load('real_Sinr.mat','*');
E = zeros(1,5);
E(1) = Sinr_real_r100 - Sinr_simul_r100.Sinr_simul;
E(2) = Sinr_real_r300 - Sinr_simul_r300.Sinr_simul;
E(3) = Sinr_real_r500 - Sinr_simul_r500.Sinr_simul;
E(4) = Sinr_real_r800 - Sinr_simul_r800.Sinr_simul;
E(5) = Sinr_real_r1200 - Sinr_simul_r1200.Sinr_simul;






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