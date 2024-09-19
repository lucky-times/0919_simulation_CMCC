function NR_plot_bs_bs_pathloss(SYS_config,UEs,activate_UE_id)
bs_bs_pathloss_sum=[];
for u_=1:length(activate_UE_id)
    bs_bs_pathloss=UEs(activate_UE_id(u_)).bs_bs_pathloss;
    bs_bs_pathloss_sum = [bs_bs_pathloss_sum bs_bs_pathloss];
end
bs_bs_pathloss_sum(isnan(bs_bs_pathloss_sum))=[];
pathloss = bs_bs_pathloss_sum;
min_ = floor(min(pathloss));
max_ = ceil(max(pathloss));
range = min_:1:max_;
h = hist(pathloss,range);
pathloss_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    plot(range,pathloss_CDF);
    xlabel('BS to BS Pathloss（dB）');
    ylabel('CDF');
end

pathloss = sort(pathloss);
index = 1:length(pathloss)/100:length(pathloss);
index = round(index);
R6 = pathloss(index);
R6 = [R6 max(pathloss)];
fid=fopen('.\calibration\R6.txt','wt');
fprintf(fid,'%f\n',R6);
fprintf(fid,'该文档记录统计到的BS到BS路损\n');
fclose(fid);

