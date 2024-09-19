clear all;
clc;
R = 6371e3; % 地球半径，单位：米
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

src_pos = [119.8395, 31.7128];
pos_1 = [119.84, 31.7136];
d = calculate_diatance(pos_1(2), pos_1(1))
x = pos - src_pos;
x = x * R;
scatter(x(:, 1), x(:, 2))
hold on;
scatter(0, 0, 'r')
hold off;


%% 
% 示例使用
% latA = ; % A点纬度
% lonA = ; % A点经度
latB = 31.7136; % B点纬度
lonB = 119.84; % B点经度

[x, y] = latlon_to_xy(latA, lonA, latB, lonB);
disp(['x: ', num2str(x), ' meters']);
disp(['y: ', num2str(y), ' meters']);


%% 
178*1e3/4/3600



