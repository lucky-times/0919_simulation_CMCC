clear all;
close all;
clc
%% 加载UAV信息
str = '7.4.2多层窄波束配置--蜂窝三扇区组网/7.4.2多层窄波束配置--蜂窝三扇区组网-1km站间距-飞行高度100米-加扰50%.csv';
data = readtable(str);
info = data(:, [3,4,5,10,8]);
info = table2array(info);
temp = info((info(:, 1)> 85 & info(:, 1) <= 135), :);% 高度、经度、纬度、SINR值

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
title('UAV h=100')
%计算UAV轨迹各个点到飞行圆心的距离
d = zeros(length(temp),1);
for i = 1:length(temp)
    d(i) = calculate_diatance(temp(i, 3), temp(i, 2));
end
hold off;

%% 绘制100m飞行高度，5种半径下的SINR值，100m,300m,500m,800m,1200m
r_100 = 105:170;
r_300 = 217:367;
r_500 = 453:675;
r_800 = 750:1108;
r_1200 = 1370:2080;


% plot_Sinr = (SINR2 - mean(SINR2))*scale_factor + mean(SINR2);
% plot(1:length(UEs), plot_Sinr,LineWidth=2);%+10+sqrt(10)*randn(size(SINR2))
% xlim([0,length(UEs)]);
% ylim(y_lim);
% grid on;
% set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
% set(gca, 'LineWidth', 1);
% ax = gca; % 获取当前坐标轴
% ax.Box = 'off'; % 关闭上边和右边的坐标轴
% 
% if(SYS_config.labelMode == 'English')
%     xlabel('UAV flight sequence number')
%     ylabel('SINR/dB')
%     str = sprintf('SINR Calculation Based on Paper Formula, h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
% elseif(SYS_config.labelMode == 'Chinese')
%     xlabel('无人机飞行序列号')
%     ylabel('SINR/(dB)');
%     str = sprintf('根据论文公式计算SINR, h = %dm, r = %dm', SYS_config.UE_height, SYS_config.UE_r);
% end  
% title(str)
% 
% str2 = sprintf('Sinr_simul_r%d.mat', SYS_config.UE_r);
labelMode = 'English'; %Chinese  English

figure
plot(1:length(r_100), temp(r_100, 4),LineWidth=2);
xlim([1,length(r_100)])
ylim([-8,15])

grid on;
set(gca, 'GridLineStyle', '--', 'GridColor', 'k', 'GridAlpha', 0.3, 'LineWidth', 0.5);
set(gca, 'LineWidth', 1);
ax = gca; % 获取当前坐标轴
ax.Box = 'off'; % 关闭上边和右边的坐标轴

if(labelMode == 'English')
    xlabel('UAV flight sequence number')
    ylabel('SINR/dB')
    str = sprintf('SINR Calculation Based on Field Data, h = 100m, r = 100m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 100m, r = 100m');
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
    str = sprintf('SINR Calculation Based on Field Data, h = 100m, r = 300m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 100m, r = 300m');
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
    str = sprintf('SINR Calculation Based on Field Data, h = 100m, r = 500m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 100m, r = 500m');
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
    str = sprintf('SINR Calculation Based on Field Data, h = 100m, r = 800m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 100m, r = 800m');
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
    str = sprintf('SINR Calculation Based on Field Data, h = 100m, r = 1200m');
elseif(labelMode == 'Chinese')
    xlabel('无人机飞行序列号')
    ylabel('SINR/(dB)');
    str = sprintf('实测数据SINR, h = 100m, r = 1200m');
end  
title(str)

save('real_Sinr.mat', "Sinr_real_r100","Sinr_real_r300","Sinr_real_r500","Sinr_real_r800","Sinr_real_r1200");
%% 计算核心橙均分
a = [2 2 2 2 2 3 2 1 1 2];
b = [99 86 89 88 97 99 88 99 90 95];
c = a*b'/sum(a)



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
