% 绘制蜂窝网络场景下，uav高度50m，距中心基站距离d=50m，theta角度0：0.3：360°分布，站间距为200m，干扰规律绘制图

close
clear all;
d = 200;
h = 50;
d_uav = 50;
%生成干扰基站位置，蜂窝网络
pos = zeros(7,2);
for i = 2:7
    pos(i, :) = [d*cosd(60*(i-2)), d*sind(60*(i-2))];
end
theta = 0:0.01:2*pi;
pos_uav = [d_uav*cos(theta); d_uav*sin(theta)];
alpha1 = atand(h/d_uav);
interfere = zeros(6,length(theta));
%计算每个干扰站对所有UAV位置点的干扰
for i = 2:7
    distance1 = sqrt((pos(i,1)-pos_uav(1,:)).^2 + (pos(i,2)-pos_uav(2,:)).^2);
    alpha2 = atand(h./distance1);
    distance2 = sqrt(pos_uav(1,:).^2 + pos_uav(2,:).^2);
    n = distance1./distance2;
    interfere(i-1, :) = 10.^((-2.2)*log10(cosd(alpha2)/cos(alpha1).*n));
end
sinr = 1./sum(interfere,1);%求出所有干扰站干扰和，并取倒数，计算SINR
plot(theta, sinr, 'r',LineStyle="-", LineWidth=2.5);
xlabel('Polar angle/rad','Interpreter','latex');
ylabel('SINR');
xlim([-0.1, 2*pi+0.1])
title('UAV h:50m, radius:50m, d between sits:200m, Macrocellulars, 6 interfering sites');
defaultAxes()
removetrTicks()
