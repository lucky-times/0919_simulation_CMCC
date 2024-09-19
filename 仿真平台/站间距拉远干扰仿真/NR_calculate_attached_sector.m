function [eNodeB_sector_id] = NR_calculate_attached_sector(SYS_config, BS_pos, UE_pos)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
m = size(BS_pos, 1);
d = zeros(m,1);
for i = 1:m
    d(i) = pdist2(BS_pos(i,:), UE_pos);
end
[~,  index] = sort(d);
eNodeB_id = 0;
for i = 1:m
    alpha = atand((SYS_config.UE_height - 25) / d(index(i)));
    if(alpha <= (SYS_config.tilt+24) && alpha > SYS_config.tilt)
        eNodeB_id = index(i);
        break;
    end
end
if eNodeB_id ==0
    error("归属基站判别出错！")
end

attached_eNodeB_pos = BS_pos(eNodeB_id,:);
angle = atan2(UE_pos(2)-attached_eNodeB_pos(2), UE_pos(1)-attached_eNodeB_pos(1))* 180 / pi;
if angle<0
    angle = angle + 360;
end
if angle == 0
    angle = 1;
end
m = ceil(angle/120);
eNodeB_sector_id =  (eNodeB_id-1)*3 + m;

end