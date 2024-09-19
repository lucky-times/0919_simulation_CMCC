% 生成信道所需参数
function [sigma_ASD,sigma_ZSD,sigma_DS,sigma_SF,sigma_ASA,sigma_ZSA,N_cluster,N_ray,...
          c_ASD,m_ZSD,Dsp,Pcs,m_ZSD_offset,c_ZSA,c_ASA,m_xpr,s_xpr]=generate_para(scene,d2D)

% 输入
% scene：场景，1为UMi，2为UMa
% d2D：基站底到移动台底的距离

% 输出
% sigma_ASD：水平离开角角度扩展
% sigma_ZSD：垂直离开角角度扩展
% sigma_DS：时延扩展
% sigma_SF:阴影衰落
% sigma_ASA：水平到达角角度扩展
% sigma_ZSA：垂直到达角角度扩展
% N_cluster：簇数
% N_ray：射线数
% c_ASD：簇的水平离开角角度扩展
% m_ZSD：垂直离开角角度扩展的mu值
% Dsp：时延缩放参数
% Pcs：每一簇的阴影衰落
% m_ZSD_offset：垂直离开角角度扩展的mu值偏移量

if scene == 1,%UMi,
    % 参数设置（都是NLOS情况）
    paraset=[-6.89,0.54,1.41,0.17,4,0,0,-0.7,0,-0.5,0.5,3,8.0,3,19,20,10,7,3,0,0,0.2,-0.4,0.4,0.5,0,0,0,0.15,1.84,0.16,0.88,22];
    
    % 垂直离开角角度扩展
    m_ZSD=max(-0.5, -2.1*(d2D/1000) +0.9);
    e_ZSD=0.6;
    m_ZSD_offset= -10^(-0.55*log10(max(10, d2D))+1.6);
else
    paraset=[-6.44,0.39,1.41,0.28,6,0.4,-0.6,-0.4,0,-0.5,0.5,2.3,7,3,20,20,2,7,3,0.4,0,1,0,0,0.6,-0.1,0,-0.4,0,0.11,1.87,0.16,1.26,15];
    % NLOS场景
    m_ZSD=max(-0.5, -2.1*(d2D/1000) +0.9);
    e_ZSD=0.49;
    m_ZSD_offset = -10^(-0.62*log10(max(10, d2D))+1.93);
end

m_ds=paraset(1);    
e_ds=paraset(2);    
m_ASD=paraset(3);    
e_ASD=paraset(4);    
s_SF=paraset(5); 
e_ASA=paraset(29);
m_ASA=paraset(30);
e_ZSA=paraset(31);
m_ZSA=paraset(32);

% 参数相关性
r_ASD_DS=paraset(6);    
r_ASD_SF=paraset(7);    
r_DS_SF=paraset(8);    
r_ZSD_SF=paraset(9);    
r_ZSD_DS=paraset(10);    
r_ZSD_ASD=paraset(11);   

r_ASA_ASD=paraset(20);
r_ASA_ZSD=paraset(21);
r_ASA_ZSA=paraset(22);
r_ASA_SF=paraset(23);
r_ASA_DS=paraset(24);
r_ZSA_ASD=paraset(25);
r_ZSA_ZSD=paraset(26);
r_ZSA_SF=paraset(27);
r_ZSA_DS=paraset(28);

% 时延缩放参数
Dsp=paraset(12); 
% 天线交叉极化比
m_xpr=paraset(13);    
s_xpr=paraset(14);    
N_cluster=paraset(15);    
N_ray=paraset(16);  
% 簇的水平离开角角度扩展
c_ASD=paraset(17);  

% 簇的垂直到达角角度扩展
c_ZSA=paraset(18); 
c_ASA=paraset(33);
% 每一簇的阴影衰落
Pcs=paraset(19);    

% 大尺度参数相关矩阵
 
A = [1           r_ZSD_ASD         r_ASA_ASD        r_ZSA_ASD         r_ASD_SF      r_ASD_DS;
     r_ZSD_ASD      1              r_ASA_ZSD        r_ZSA_ZSD         r_ZSD_SF      r_ZSD_DS;
     r_ASA_ASD   r_ASA_ZSD             1            r_ASA_ZSA         r_ASA_SF      r_ASA_DS;
     r_ZSA_ASD   r_ZSA_ZSD         r_ASA_ZSA            1             r_ZSA_SF      r_ZSA_DS;
     r_ASD_SF    r_ZSD_SF          r_ASA_SF         r_ZSA_SF             1          r_DS_SF;
     r_ASD_DS    r_ZSD_DS          r_ASA_DS         r_ZSA_DS          r_DS_SF          1;];
     
 
w0 = randn(6,1);
x  = A^(0.5)*w0;

sigma_ASD = 10^(e_ASD*x(1)+m_ASD);
sigma_ZSD = 10^(e_ZSD*x(2)+m_ZSD);
sigma_DS = 10^(e_ds*x(6)+m_ds);
sigma_SF = 10^((s_SF*x(5)/10)); 
sigma_ASA = 10^(e_ASA*x(3)+m_ASA);
sigma_ZSA = 10^(e_ZSA*x(4)+m_ZSA);

end

