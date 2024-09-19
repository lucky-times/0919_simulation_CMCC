function pos = latlon_to_xy(B)
%UNTITLED 计算以地球坐标系中latA, lonA为原点，B相对于A的x,y坐标，球面坐标系转化为直角坐标系
%   输入：B，代表待转换点的经纬度
%   输出：pos包含x,y
%   x, B相对于A的直角坐标系的横坐标，y, B相对于A的直角坐标系的横坐标， 单位：m
latA = 31.7128;
lonA = 119.8395;
% 将经纬度转换为弧度
latA = deg2rad(latA);
lonA = deg2rad(lonA);
latB = deg2rad(B(2));
lonB = deg2rad(B(1));

% 地球半径
R = 6371e3; % 单位：米

% 计算B点相对于A点的球面坐标
dlat = latB - latA;
dlon = lonB - lonA;

% 平面近似公式
x = dlon * cos((latA + latB) / 2) * R;
y = dlat * R;
pos = [x ,y];
end
