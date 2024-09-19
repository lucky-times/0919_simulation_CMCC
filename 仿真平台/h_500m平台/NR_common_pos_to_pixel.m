function [ pos_pixel,pos_pixel_exact] = NR_common_pos_to_pixel( pos,       roi_min, data_res)
% position，将绝对值位置转换到像素位置，和另一个对应
% 将像素点位置转换到位置绝对值
% 输入参数：
% pos：实际的位置
% roi_min：ROI最左侧的和最下侧值
% data_res：数据分辨率
%
% 输出参数：
% pos_pixel：像素点的位置


pos_pixel(:,1) = floor((pos(:,1)-roi_min(1))/data_res)+1;
pos_pixel(:,2) = floor((pos(:,2)-roi_min(2))/data_res)+1;

pos_pixel_exact(:,1) = ((pos(:,1)-roi_min(1))/data_res)+1;
pos_pixel_exact(:,2) = ((pos(:,2)-roi_min(2))/data_res)+1;

