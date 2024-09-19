%function NR_dl_ul_gain(SYS_config,UEs,activate_UE_id)
function NR_user_PL_gain(SYS_config,user_macroscopic_pathloss,dl_user_macroscopic_pathloss)
%NR_UL_DL_ 此处显示有关此函数的摘要
%  保存
% dl_ul_gain_tx_sum=[];
% for u_=1:length(activate_UE_id)
%     dl_ul_gain_tx=UEs(activate_UE_id(u_)).dl_ul_tx;
%     dl_ul_gain_tx_sum = [ dl_ul_gain_tx_sum dl_ul_gain_tx];
% end
% dl_ul_gain_tx_sum (isnan(dl_ul_gain_tx_sum ))=[];
% dl_ul_tx_gain = dl_ul_gain_tx_sum ;
min_ = floor(min(user_macroscopic_pathloss));
max_ = ceil(max(user_macroscopic_pathloss));
range = min_:1:max_;
h = hist(user_macroscopic_pathloss,range);
user_macroscopic_pathloss_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    subplot(1,2,1);
    plot(range,user_macroscopic_pathloss_CDF);
    hold on;
    xlabel('user_macroscopic_pathloss');
    ylabel('CDF');
end

user_macroscopic_pathloss = sort(user_macroscopic_pathloss);  %-20
index = 1:length(user_macroscopic_pathloss)/100:length(user_macroscopic_pathloss);
index = round(index);
R = user_macroscopic_pathloss(index);
R = [R max(user_macroscopic_pathloss)];
fid=fopen('.\calibration\user_macroscopic_pathloss.txt','wt');
fprintf(fid,'%f\n',R);
fclose(fid);

min_1 = floor(min(dl_user_macroscopic_pathloss));
max_1 = ceil(max(dl_user_macroscopic_pathloss));
range2 = min_1:1:max_1;
h2 = hist(dl_user_macroscopic_pathloss,range2);
dl_user_macroscopic_pathloss_CDF = cumsum(h2)/sum(h2);
if SYS_config.photo
    figure;
    subplot(1,2,2);
    plot(range2,dl_user_macroscopic_pathloss_CDF);
    hold on;
    xlabel('dl_user_macroscopic_pathloss');
    ylabel('CDF');
end

dl_user_macroscopic_pathloss = sort(dl_user_macroscopic_pathloss);
index2 = 1:length(dl_user_macroscopic_pathloss)/100:length(dl_user_macroscopic_pathloss);
index2 = round(index2);
R2 = dl_user_macroscopic_pathloss(index2);
R2 = [R2 max(dl_user_macroscopic_pathloss)];
fid=fopen('.\calibration\dl_user_macroscopic_pathloss.txt','wt');
fprintf(fid,'%f\n',R2);
fclose(fid);


end

