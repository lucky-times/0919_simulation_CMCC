%����:   scene: ����ѡ�� 1 Ϊ UMi, 2 Ϊ UMa
%          d2D: С���뾶
%            U: �ƶ���������
%            S: ��վ��������
%���:       H:  �ŵ�����
function H = channel_3DMIMO(scene,d2D,U,S)

clusterE = pi/2+(2*rand-1)*45/180*pi; %  phi_zod
clusterD = pi/2+(2*rand-1)*60/180*pi; %  phi_aod

[sigma_ASD,sigma_ZSD,sigma_DS,sigma_SF,sigma_ASA,sigma_ZSA,N_cluster,N_ray,...
 c_ASD,m_ZSD,Dsp,Pcs,m_ZSD_offset,c_ZSA,c_ASA,m_xpr,s_xpr]=generate_para(scene,d2D);

[pha_aod_n_m,pha_zod_n_m,pha_aoa_n_m,pha_zoa_n_m,p]=angles(clusterE,clusterD,...
 sigma_ASD,sigma_ZSD,sigma_ASA,sigma_ZSA,sigma_DS,N_cluster,N_ray,c_ASD,c_ASA,c_ZSA,m_ZSD,Dsp,Pcs,m_ZSD_offset);

% �������߽��漫����
m_XPR = m_xpr;
s_XPR = s_xpr;
XPR = 10.^(normrnd(m_XPR,s_XPR,U,S)/10);
alpha = pi/4;

% �������ּ�����λ
phi_VV = pi.*rand(N_ray,N_cluster);
phi_VH = pi.*rand(N_ray,N_cluster);
phi_HV = pi.*rand(N_ray,N_cluster);
phi_HH = pi.*rand(N_ray,N_cluster);

% �����ŵ�ϵ��
t = [1:100];
CarrierFrequency = 2e9;
wavelength = 3e8/CarrierFrequency;
ds = 0.5*wavelength;
du = 0.5*wavelength;
k_CONST = 2*pi/wavelength;       % ��ʸ
v = 10;                          % m/s
theta_v = pi*rand;
% �ŵ������ʼ��
H = zeros(U,S,length(t));
h_m(t) = zeros(size(t));
h_n(t) = zeros(size(t));
h(t) = zeros(size(t));
for u=1:U
    for s=1:S
        for n=1:N_cluster
            for m=1:N_ray
            h_m(t) = sqrt(p(n)/N_ray)*[cos(alpha)*sin(pha_zoa_n_m(m,n))+sin(alpha)*sin(pha_zoa_n_m(m,n))*cos(pha_aoa_n_m(m,n)),cos(alpha)*cos(pha_zoa_n_m(m,n))]* ...
                     [exp(1i*phi_VV(m,n)),sqrt(XPR(u,s))*exp(1i*phi_VH(m,n));sqrt(XPR(u,s))*exp(1i*phi_HV(m,n)),exp(1i*phi_HH(m,n))]*...
                     [cos(alpha)*sin(pha_zod_n_m(m,n))+sin(alpha)*sin(pha_zod_n_m(m,n))*cos(pha_aod_n_m(m,n));cos(alpha)*cos(pha_zod_n_m(m,n))]*...
                     exp(1i*(k_CONST*ds*sin(pha_aod_n_m(m,n))*sin(pha_zod_n_m(m,n))))*...
                     exp(1i*(k_CONST*du*sin(pha_aoa_n_m(m,n))*sin(pha_zoa_n_m(m,n)))).*...
                     exp(1i*k_CONST)*v*sin(pha_zoa_n_m(m,n)*cos(pha_aoa_n_m(m,n)-theta_v).*t);
            h_n(t) = h_n(t)+h_m(t);
            end
            h(t) = h(t)+h_n(t);
        end
        H(u,s,:) = h(t);
    end
end
end
        
