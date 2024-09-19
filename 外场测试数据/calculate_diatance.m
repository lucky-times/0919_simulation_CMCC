function [distance] = calculate_diatance(point1_lat,point1_lon)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
% 将经纬度转换为弧度
point2_lat = 31.7128;
point2_lon = 119.84;
point1_lat_rad = deg2rad(point1_lat);
point1_lon_rad = deg2rad(point1_lon);
point2_lat_rad = deg2rad(point2_lat);
point2_lon_rad = deg2rad(point2_lon);

% 计算两点之间的距离
R = 6371e3; % 地球半径，单位：千米
dlat = point2_lat_rad - point1_lat_rad;
dlon = point2_lon_rad - point1_lon_rad;
a = sin(dlat/2)^2 + cos(point1_lat_rad) * cos(point2_lat_rad) * sin(dlon/2)^2;
c = 2 * atan2(sqrt(a), sqrt(1-a));
distance = R * c;

% 输出结果
% fprintf('两点之间的距离为：%.2f km\n', distance);
end