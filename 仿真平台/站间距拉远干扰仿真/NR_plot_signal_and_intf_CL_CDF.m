function NR_plot_signal_and_intf_CL_CDF(signal_CL,interfering_CL,another_interfering_CL,activate_UE_id)
% 保存每个系统CL的cdf以及服务信号的CL
% 输入参数：
% signal_CL：服务信号的耦合损耗
% interfering_CL：同频干扰信号的CL
% another_interfering_CL：邻频干扰信号的CL
% UEs：UE列表

n_UE = length(activate_UE_id);
% 确定待保存服务信号CL的索引（为了获得101个数据便于校准）
s_CL = -signal_CL';
s_CL = sort(s_CL);
index = 1:length(s_CL)/100:length(s_CL);
index = round(index);
% 保存服务信号
R_s_CL = s_CL(index);
R_s_CL = [R_s_CL max(s_CL)];
fid=fopen('.\calibration\R_s_CL.txt','wt');
fprintf(fid,'%f\n',R_s_CL);
fclose(fid);
% 确定待保存同频干扰信号CL的索引
CL = reshape(interfering_CL,1,n_UE*size(interfering_CL,2));
CL = -CL;
CL = sort(CL);
index = 1:length(CL)/100:length(CL);
index = round(index);
% 保存同频的CL
R_CL = CL(index);
R_CL = [R_CL max(CL)];
fid=fopen('.\calibration\R_CL.txt','wt');
fprintf(fid,'%f\n',R_CL);
fclose(fid);
% 确定待保存邻频干扰信号CL的索引
another_CL = reshape(another_interfering_CL,1,n_UE*size(another_interfering_CL,2));
another_CL = -another_CL;
another_CL = sort(another_CL);
index = 1:length(another_CL)/100:length(another_CL);
index = round(index);
% 保存邻频的CL
R_another_CL = another_CL(index);
R_another_CL = [R_another_CL max(another_CL)];
fid=fopen('.\calibration\R_another_CL.txt','wt');
fprintf(fid,'%f\n',R_another_CL);
fclose(fid);
end