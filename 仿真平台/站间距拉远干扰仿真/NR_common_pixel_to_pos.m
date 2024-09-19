function [ pos ]                      = NR_common_pixel_to_pos( pos_pixel, roi_min, data_res)
% 将像素点位置转换到位置绝对值
% 输入参数：
% pos_pixel：像素点的位置
% roi_min：ROI最左侧的和最下侧值
% data_res：数据分辨率
%
% 输出参数：
% pos：实际的位置

pos(:,1) = (pos_pixel(:,1)-1)*data_res+roi_min(1);
pos(:,2) = (pos_pixel(:,2)-1)*data_res+roi_min(2);
