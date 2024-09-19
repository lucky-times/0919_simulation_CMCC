% ����AoD,ZoD,AoA,ZoA
function [pha_aod_n_m,pha_zod_n_m,pha_aoa_n_m,pha_zoa_n_m,p]=angles(clusterE,clusterD,...
          sigma_ASD,sigma_ZSD,sigma_ASA,sigma_ZSA,sigma_DS,N_cluster,N_ray,c_ASD,c_ASA,c_ZSA,m_ZSD,Dsp,Pcs,m_ZSD_offset)
% ����
% sigma_ASD��ˮƽ�뿪�ǽǶ���չ
% sigma_ZSD����ֱ�뿪�ǽǶ���չ
% sigma_DS��ʱ����չ
% N_cluster������
% N_ray���Ӿ���
% c_ASD���ص�ˮƽ�뿪�ǽǶ���չ
% m_ZSD����ֱ�뿪�ǽǶ���չ��muֵ
% Dsp��ʱ�����Ų���
% Pcs��ÿһ�ص���Ӱ˥��
% m_ZSD_offset����ֱ�뿪�ǽǶ���չ��muֵƫ����

% ���
% pha_aod_n_m��ÿһ�Ӿ���ˮƽ�뿪��
% pha_zod_n_m��ÿһ�Ӿ��Ĵ�ֱ�뿪��
% p��N���ྶ�����Ĺ�һ������

% ˮƽ�뿪�Ǻʹ�ֱ�뿪��
if N_cluster == 20,
    cnst_AoD = 1.289;  % ��������C
    cnst_ZoD = 1.178;
    cnst_AoA = 1.289;
    cnst_ZoA = 1.178;
else                   % N_cluster=12
    cnst_AoD = 1.146;
    cnst_ZoD = 1.104;
    cnst_AoA = 1.146;
    cnst_ZoA = 1.104;
end

% �ص�ˮƽ�뿪�Ǻʹ�ֱ�뿪��
c_AoD = c_ASD;
c_ZoD = 3*10^(m_ZSD)/8;
c_AoA = c_ASA;
c_ZoA = c_ZSA;

% ȷ��N���ྶ���������ʱ��
x=rand(1,N_cluster);
tau1 = -Dsp*sigma_DS*log(x);
tau = sort(tau1,'descend');

% ȷ��N���ྶ���������ƽ������
z = Pcs*randn(1,N_cluster);
p1 = exp(-tau*(Dsp-1)/(Dsp*sigma_DS)).*10.^(-z/10);
p = p1/sum(p1); % ��һ����ʹ����֮��Ϊ1

MaxP= max(p);
% ���� AoD  p.26
pha_aod1 = 2/1.4/cnst_AoD*sigma_ASD*sqrt(-log(p/MaxP));

Xn = sign(randn(1,N_cluster)); % ֵΪ1��-1�ľ��ȷֲ�
Yn = sigma_ASD/7*randn(1,N_cluster);  % Yn�ֲ���N(0,sigma_ASD^2/7^2)
pha_aod =  Xn.*pha_aod1+Yn;
pha_aod_n_m = (repmat(pha_aod,N_ray,1)+c_AoD*randn(N_ray,N_cluster))/180*pi+clusterD;
pha_aod_n_m = real(pha_aod_n_m);

pha_aod_n_m=mod(pha_aod_n_m,2*pi);
% ʹAoD��ΧΪ0��Pi
index_a = find(pha_aod_n_m>pi);
pha_aod_n_m(index_a)=2*pi-pha_aod_n_m(index_a);

% ���� ZoD
pha_zod1 = 1/cnst_ZoD*sigma_ZSD*(-log(p/MaxP));

Xn = sign(randn(1,N_cluster));
Yn = sigma_ZSD/7*randn(1,N_cluster);
pha_zod =  Xn.*pha_zod1+Yn;
pha_zod_n_m = (repmat(pha_zod,N_ray,1)+c_ZoD*randn(N_ray,N_cluster)+m_ZSD_offset)/180*pi+clusterE;
pha_zod_n_m = real(pha_zod_n_m);

pha_zod_n_m=mod(pha_zod_n_m,2*pi);
index_z = find(pha_zod_n_m>pi);
pha_zod_n_m(index_z)=2*pi-pha_zod_n_m(index_z);

% ���� AoA
pha_aoa1 = 2/1.4/cnst_AoA*sigma_ASA*sqrt(-log(p/MaxP));

Xn = sign(randn(1,N_cluster)); 
Yn = sigma_ASA/7*randn(1,N_cluster);  
pha_aoa =  Xn.*pha_aoa1+Yn;
pha_aoa_n_m = (repmat(pha_aoa,N_ray,1)+c_AoA*randn(N_ray,N_cluster))/180*pi+clusterD;
pha_aoa_n_m = real(pha_aoa_n_m);

pha_aoa_n_m=mod(pha_aoa_n_m,2*pi);
% ʹAoA��ΧΪ0��Pi
index_aoa = find(pha_aoa_n_m>pi);
pha_aoa_n_m(index_aoa)=2*pi-pha_aoa_n_m(index_aoa);

% ���� ZoA
pha_zoa1 = 1/cnst_ZoA*sigma_ZSA*(-log(p/MaxP));

Xn = sign(randn(1,N_cluster));
Yn = sigma_ZSA/7*randn(1,N_cluster);
pha_zoa =  Xn.*pha_zoa1+Yn;
pha_zoa_n_m = (repmat(pha_zoa,N_ray,1)+c_ZoA*randn(N_ray,N_cluster))/180*pi+clusterE;
pha_zoa_n_m = real(pha_zoa_n_m);

pha_zoa_n_m=mod(pha_zoa_n_m,2*pi);
index_zoa = find(pha_zoa_n_m>pi);
pha_zoa_n_m(index_zoa)=2*pi-pha_zoa_n_m(index_zoa);
