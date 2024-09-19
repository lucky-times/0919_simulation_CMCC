% 生成AoD,ZoD,AoA,ZoA
function [pha_aod_n_m,pha_zod_n_m,pha_aoa_n_m,pha_zoa_n_m,p]=angles(clusterE,clusterD,...
          sigma_ASD,sigma_ZSD,sigma_ASA,sigma_ZSA,sigma_DS,N_cluster,N_ray,c_ASD,c_ASA,c_ZSA,m_ZSD,Dsp,Pcs,m_ZSD_offset)
% 输入
% sigma_ASD：水平离开角角度扩展
% sigma_ZSD：垂直离开角角度扩展
% sigma_DS：时延扩展
% N_cluster：簇数
% N_ray：子径数
% c_ASD：簇的水平离开角角度扩展
% m_ZSD：垂直离开角角度扩展的mu值
% Dsp：时延缩放参数
% Pcs：每一簇的阴影衰落
% m_ZSD_offset：垂直离开角角度扩展的mu值偏移量

% 输出
% pha_aod_n_m：每一子径的水平离开角
% pha_zod_n_m：每一子径的垂直离开角
% p：N条多径分量的归一化功率

% 水平离开角和垂直离开角
if N_cluster == 20,
    cnst_AoD = 1.289;  % 缩放因子C
    cnst_ZoD = 1.178;
    cnst_AoA = 1.289;
    cnst_ZoA = 1.178;
else                   % N_cluster=12
    cnst_AoD = 1.146;
    cnst_ZoD = 1.104;
    cnst_AoA = 1.146;
    cnst_ZoA = 1.104;
end

% 簇的水平离开角和垂直离开角
c_AoD = c_ASD;
c_ZoD = 3*10^(m_ZSD)/8;
c_AoA = c_ASA;
c_ZoA = c_ZSA;

% 确定N条多径分量的随机时延
x=rand(1,N_cluster);
tau1 = -Dsp*sigma_DS*log(x);
tau = sort(tau1,'descend');

% 确定N条多径分量的随机平均功率
z = Pcs*randn(1,N_cluster);
p1 = exp(-tau*(Dsp-1)/(Dsp*sigma_DS)).*10.^(-z/10);
p = p1/sum(p1); % 归一化，使功率之和为1

MaxP= max(p);
% 生成 AoD  p.26
pha_aod1 = 2/1.4/cnst_AoD*sigma_ASD*sqrt(-log(p/MaxP));

Xn = sign(randn(1,N_cluster)); % 值为1和-1的均匀分布
Yn = sigma_ASD/7*randn(1,N_cluster);  % Yn分布：N(0,sigma_ASD^2/7^2)
pha_aod =  Xn.*pha_aod1+Yn;
pha_aod_n_m = (repmat(pha_aod,N_ray,1)+c_AoD*randn(N_ray,N_cluster))/180*pi+clusterD;
pha_aod_n_m = real(pha_aod_n_m);

pha_aod_n_m=mod(pha_aod_n_m,2*pi);
% 使AoD范围为0到Pi
index_a = find(pha_aod_n_m>pi);
pha_aod_n_m(index_a)=2*pi-pha_aod_n_m(index_a);

% 生成 ZoD
pha_zod1 = 1/cnst_ZoD*sigma_ZSD*(-log(p/MaxP));

Xn = sign(randn(1,N_cluster));
Yn = sigma_ZSD/7*randn(1,N_cluster);
pha_zod =  Xn.*pha_zod1+Yn;
pha_zod_n_m = (repmat(pha_zod,N_ray,1)+c_ZoD*randn(N_ray,N_cluster)+m_ZSD_offset)/180*pi+clusterE;
pha_zod_n_m = real(pha_zod_n_m);

pha_zod_n_m=mod(pha_zod_n_m,2*pi);
index_z = find(pha_zod_n_m>pi);
pha_zod_n_m(index_z)=2*pi-pha_zod_n_m(index_z);

% 生成 AoA
pha_aoa1 = 2/1.4/cnst_AoA*sigma_ASA*sqrt(-log(p/MaxP));

Xn = sign(randn(1,N_cluster)); 
Yn = sigma_ASA/7*randn(1,N_cluster);  
pha_aoa =  Xn.*pha_aoa1+Yn;
pha_aoa_n_m = (repmat(pha_aoa,N_ray,1)+c_AoA*randn(N_ray,N_cluster))/180*pi+clusterD;
pha_aoa_n_m = real(pha_aoa_n_m);

pha_aoa_n_m=mod(pha_aoa_n_m,2*pi);
% 使AoA范围为0到Pi
index_aoa = find(pha_aoa_n_m>pi);
pha_aoa_n_m(index_aoa)=2*pi-pha_aoa_n_m(index_aoa);

% 生成 ZoA
pha_zoa1 = 1/cnst_ZoA*sigma_ZSA*(-log(p/MaxP));

Xn = sign(randn(1,N_cluster));
Yn = sigma_ZSA/7*randn(1,N_cluster);
pha_zoa =  Xn.*pha_zoa1+Yn;
pha_zoa_n_m = (repmat(pha_zoa,N_ray,1)+c_ZoA*randn(N_ray,N_cluster))/180*pi+clusterE;
pha_zoa_n_m = real(pha_zoa_n_m);

pha_zoa_n_m=mod(pha_zoa_n_m,2*pi);
index_zoa = find(pha_zoa_n_m>pi);
pha_zoa_n_m(index_zoa)=2*pi-pha_zoa_n_m(index_zoa);
