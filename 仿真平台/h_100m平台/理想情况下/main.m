% 只考虑基站的干扰，没有考虑三扇区，利用ni和alpha计算出SINR


clear all
ISD = 300;
n_rings = 2;
eNodeB_pos = creat_eNodeB(ISD, n_rings);
h = 50;
r = 88; %极径
step = 0.1;
theta = 0:step:360-step;%极角
UE_positions = [r*cosd(theta);r*sind(theta)];

d = zeros(size(eNodeB_pos, 2), length(theta));
alpha = atand(h/r);
for i = 1:size(eNodeB_pos, 2)
    for j = 1:length(theta)
        d(i, j) = sqrt((eNodeB_pos(1,i)-UE_positions(1,j))^2 + (eNodeB_pos(2,i)-UE_positions(2,j))^2);
    end
end

n = zeros(size(d));
for i = 1:size(eNodeB_pos, 2)
    for j = 1:length(theta)
        n(i, j) = d(i,j)/d(11, j);
    end
end
% 删除主服务小区的距离比
n(11, :) = [];

inter = zeros(size(n));
for i = 1:size(inter, 1)
    for j = 1:size(inter, 2)
        inter(i, j) = 10^(-2.2*log10(cosd(alpha)/cosd(alpha/n(i,j))*n(i,j)));
    end
end
inter_sum = sum(inter, 1);
SINR = 1./inter_sum;
plot(theta, SINR, LineWidth=2)
xlabel('极角')
ylabel('SINR')
xlim([-5, 365])



