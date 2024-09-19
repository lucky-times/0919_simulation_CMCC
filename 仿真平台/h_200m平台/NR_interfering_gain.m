%function NR_dl_ul_gain(SYS_config,UEs,activate_UE_id)
function NR_interfering_gain(SYS_config,UEs,activate_UE_id)
%NR_UL_DL_ 此处显示有关此函数的摘要
%  保存
dl_ul_UE_gain_sum=[];
dl_ul_BS_gain_sum=[];
for u_=1:length(activate_UE_id)
    dl_ul_gain_UE=UEs(activate_UE_id(u_)).dl_ul_UE_gain;
    dl_ul_gain_BS=UEs(activate_UE_id(u_)).dl_ul_BS_gain;
    dl_ul_UE_gain_sum = [ dl_ul_UE_gain_sum dl_ul_gain_UE];
    dl_ul_BS_gain_sum = [ dl_ul_BS_gain_sum dl_ul_gain_BS];
end
dl_ul_UE_gain_sum (isnan(dl_ul_UE_gain_sum ))=[];
dl_ul_BS_gain_sum (isnan(dl_ul_BS_gain_sum ))=[];
dl_ul_UE_gain = -dl_ul_UE_gain_sum ;
dl_ul_BS_gain = dl_ul_BS_gain_sum ;
min_1 = floor(min(dl_ul_UE_gain));
max_1 = ceil(max(dl_ul_UE_gain));
range = min_1:1:max_1;
h = hist(dl_ul_UE_gain,range);
dl_ul_UE_gain_CDF = cumsum(h)/sum(h);

min_2 = floor(min(dl_ul_BS_gain));
max_2 = ceil(max(dl_ul_BS_gain));
range2 = min_2:1:max_2;
h2 = hist(dl_ul_BS_gain,range2);
dl_ul_BS_gain_CDF = cumsum(h2)/sum(h2);

if SYS_config.photo
    figure;
    subplot(2,2,1);
    plot(range,dl_ul_UE_gain_CDF);
    grid on;
    xlabel('dl ul UE gain');
    ylabel('CDF');

    subplot(2,2,2);
    plot(range2,dl_ul_BS_gain_CDF);
    grid on;
    xlabel('dl to ul BS gain');
    ylabel('CDF');
end

dl_ul_UE_gain= sort(dl_ul_UE_gain);
index = 1:length(dl_ul_UE_gain)/100:length(dl_ul_UE_gain);
index = round(index);
R1 = dl_ul_UE_gain(index);
R1 = [R1 max(dl_ul_UE_gain)];
fid=fopen('.\calibration\dl_ul_UE_gain.txt','wt');
fprintf(fid,'%f\n',R1);
fclose(fid);

dl_ul_BS_gain= sort(dl_ul_BS_gain);
index2 = 1:length(dl_ul_BS_gain)/100:length(dl_ul_BS_gain);
index2 = round(index2);
R2 = dl_ul_BS_gain(index2);
R2 = [R2 max(dl_ul_BS_gain)];
fid=fopen('.\calibration\dl_ul_BS_gain.txt','wt');
fprintf(fid,'%f\n',R2);
fclose(fid);

ul_dl_UE_gain_sum=[];
ul_dl_BS_gain_sum=[];
for u_=1:length(activate_UE_id)
    ul_dl_gain_UE=UEs(activate_UE_id(u_)).ul_dl_UE_gain';
    ul_dl_gain_BS=UEs(activate_UE_id(u_)).ul_dl_BS_gain';
    ul_dl_UE_gain_sum = [ ul_dl_UE_gain_sum ul_dl_gain_UE];
    ul_dl_BS_gain_sum = [ ul_dl_BS_gain_sum ul_dl_gain_BS];
end
ul_dl_UE_gain_sum (isnan(ul_dl_UE_gain_sum ))=[];
ul_dl_BS_gain_sum (isnan(ul_dl_BS_gain_sum ))=[];
ul_dl_UE_gain = ul_dl_UE_gain_sum ;
ul_dl_BS_gain = ul_dl_BS_gain_sum ;
min_3 = floor(min(ul_dl_UE_gain));
max_3 = ceil(max(ul_dl_UE_gain));
range3 = min_3:1:max_3;
h3 = hist(ul_dl_UE_gain,range3);
ul_dl_UE_gain_CDF = cumsum(h3)/sum(h3);

min_4 = floor(min(ul_dl_BS_gain));
max_4 = ceil(max(ul_dl_BS_gain));
range4 = min_4:1:max_4;
h4 = hist(ul_dl_BS_gain,range4);
ul_dl_BS_gain_CDF = cumsum(h4)/sum(h4);

if SYS_config.photo
    subplot(2,2,3);
    plot(range3,ul_dl_UE_gain_CDF);
    grid on;
    xlabel('ul dl UE gain');
    ylabel('CDF');

    subplot(2,2,4);
    plot(range4,ul_dl_BS_gain_CDF);
    grid on;
    xlabel('ul dl BS gain');
    ylabel('CDF');
    hold off
end

ul_dl_UE_gain= sort(ul_dl_UE_gain);
index3 = 1:length(ul_dl_UE_gain)/100:length(ul_dl_UE_gain);
index3 = round(index3);
R3 = ul_dl_UE_gain(index3);
R3 = [R3 max(ul_dl_UE_gain)];
fid=fopen('.\calibration\ul_dl_UE_gain.txt','wt');
fprintf(fid,'%f\n',R3);
fclose(fid);

ul_dl_BS_gain= sort(ul_dl_BS_gain);
index4 = 1:length(ul_dl_BS_gain)/100:length(ul_dl_BS_gain);
index4 = round(index4);
R4 = ul_dl_BS_gain(index4);
R4 = [R4 max(ul_dl_BS_gain)];
fid=fopen('.\calibration\ul_dl_BS_gain.txt','wt');
fprintf(fid,'%f\n',R4);
fclose(fid);

%被干扰链路的
dl_ul_in_UE_gain_sum=[];
dl_ul_in_BS_gain_sum=[];
for u_=1:length(activate_UE_id)
    dl_ul_in_gain_UE=UEs(activate_UE_id(u_)).interfering_UE;
    dl_ul_in_gain_BS=UEs(activate_UE_id(u_)).interfering_BS;
    dl_ul_in_UE_gain_sum = [ dl_ul_in_UE_gain_sum dl_ul_in_gain_UE];
    dl_ul_in_BS_gain_sum = [ dl_ul_in_BS_gain_sum dl_ul_in_gain_BS];
end
dl_ul_in_UE_gain_sum (isnan(dl_ul_in_UE_gain_sum ))=[];
dl_ul_in_BS_gain_sum (isnan(dl_ul_in_BS_gain_sum ))=[];
dl_ul_in_UE_gain = dl_ul_in_UE_gain_sum ;
dl_ul_in_BS_gain = dl_ul_in_BS_gain_sum ;
min_5 = floor(min(dl_ul_in_UE_gain));
max_5 = ceil(max(dl_ul_in_UE_gain));
range5 = min_5:1:max_5;
h5 = hist(dl_ul_in_UE_gain,range5);
dl_ul_in_UE_gain_CDF = cumsum(h5)/sum(h5);

min_6 = floor(min(dl_ul_in_BS_gain));
max_6 = ceil(max(dl_ul_in_BS_gain));
range6 = min_6:1:max_6;
h6 = hist(dl_ul_in_BS_gain,range6);
dl_ul_in_BS_gain_CDF = cumsum(h6)/sum(h6);

% if SYS_config.photo
%     figure();
%     subplot(2,2,1);
%     plot(range5,dl_ul_in_UE_gain_CDF);
%     grid on;
%     xlabel('dl to ul in UE gain');
%     ylabel('CDF');
% 
%     subplot(2,2,2);
%     %plot(range6,dl_ul_in_BS_gain_CDF);
%     %grid on;
%     %xlabel('dl to ul in BS gain');
%     %ylabel('CDF');
% end
% 
dl_ul_in_UE_gain= sort(dl_ul_in_UE_gain);
index5 = 1:length(dl_ul_in_UE_gain)/100:length(dl_ul_in_UE_gain);
index5 = round(index5);
R5 = dl_ul_in_UE_gain(index5);
R5 = [R5 max(dl_ul_in_UE_gain)];
fid=fopen('.\calibration\dl_ul_in_UE_gain.txt','wt');
fprintf(fid,'%f\n',R5);
fclose(fid);

dl_ul_in_BS_gain= sort(dl_ul_in_BS_gain);
index6 = 1:length(dl_ul_in_BS_gain)/100:length(dl_ul_in_BS_gain);
index6 = round(index6);
R6 = dl_ul_in_BS_gain(index6);
R6 = [R6 max(dl_ul_in_BS_gain)];
fid=fopen('.\calibration\dl_ul_in_BS_gain.txt','wt');
fprintf(fid,'%f\n',R6);
fclose(fid);

end


