function [range,output] = NR_plot_CDF(input)
% 对要求画CDF的数据进行预处理
% 输入参数
% input：待处理数据
% 输出参数
% range：CDF图中的X轴
% output：CDF图中的Y轴

min_ = floor(min(input));
max_ = ceil(max(input));
range = min_:1:max_;
h = hist(input,range);
output = cumsum(h)/sum(h);
end