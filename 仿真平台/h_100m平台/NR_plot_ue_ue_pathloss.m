function NR_plot_ue_ue_pathloss(SYS_config,UEs,activate_UE_id)
ue_ue_pathloss_sum=[];
for u_=1:length(activate_UE_id)
    ue_ue_pathloss=UEs(activate_UE_id(u_)).ue_ue_pathloss;
    ue_ue_pathloss_sum = [ue_ue_pathloss_sum ue_ue_pathloss];
end
ue_ue_pathloss_sum(isnan(ue_ue_pathloss_sum))=[];
pathloss = ue_ue_pathloss_sum;
min_ = floor(min(pathloss));
max_ = ceil(max(pathloss));
range = min_:1:max_;
h = hist(pathloss,range);
pathloss_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    plot(range,pathloss_CDF);
    xlabel('UE to UE Pathloss（dB）');
    ylabel('CDF');
end
pathloss = reshape(pathloss,1,[]);
pathloss = sort(pathloss);
index = 1:length(pathloss)/100:length(pathloss);
index = round(index);
R9 = pathloss(index);
pathloss_max = max(pathloss);
R9 = [R9 pathloss_max(1)];
fid=fopen('.\calibration\R9.txt','wt');
fprintf(fid,'%f\n',R9);
fprintf(fid,'该文档记录统计到的UE到UE路损\n');
fclose(fid);

