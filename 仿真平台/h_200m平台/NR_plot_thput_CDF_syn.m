function NR_plot_thput_CDF_syn(SYS_config,UEs,activate_UE_id)
% 绘制吞吐量cdf
% 输入参数：
% SYS_config：LTE参数系统
% UEs：UE列表
%
%一些初始化
all_SINR = [UEs.ul_dl_wideband_SINR];
ul_all_SINR = [UEs.dl_ul_wideband_SINR];
n_UE_served_per_BS = SYS_config.n_UE_served_per_BS;
n_UE = length(activate_UE_id);
bandwidth = SYS_config.bandwidth/n_UE_served_per_BS;

n_snapshot = SYS_config.UE_per_eNodeB/n_UE_served_per_BS;

% 下行
UE_SINR = all_SINR(activate_UE_id);
for i=1:length(UE_SINR)
    if UE_SINR(i)<-10
        UE_thput(i) = 0;
    elseif UE_SINR(i)>30
        UE_thput(i) = 0.6 * bandwidth/1e6*log2(1+10^(0.1*30));
    else
        UE_thput(i) = 0.6 * bandwidth/1e6*log2(1+10^(0.1*UE_SINR(i)));
    end
end

thput_per_UE = sum(UE_thput)/n_UE;
thput_per_cell = thput_per_UE*n_UE_served_per_BS;

if n_snapshot~=1
    % 例子：第一次快照基站1服务UE1和2，基站2服务UE11和12……，然后统计该快照下边缘UE的吞吐量，最后与所有快照下边缘UE的吞吐量进行平均
    UE_thput_divi = zeros(n_snapshot,n_UE/n_snapshot);
    for i = 1:n_snapshot
        index_UE = [];
        for a_ = 1:n_UE_served_per_BS
            if isempty(index_UE)
                index_UE = (i*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE;
            else
                index_UE = [index_UE,(i*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE ];
            end
        end
        UE_thput_divi(i,:) = UE_thput(index_UE);
    end
    UE_thput_divi = sort(UE_thput_divi,2);
    UE_thput_ave = sum(UE_thput_divi)/n_snapshot;
    
    index_5 = ceil(n_UE/n_snapshot*0.05);
    thput_per_UE_5 = UE_thput_ave(index_5);
else
    UE_thput_sort = sort(UE_thput);
    index_5 = ceil(n_UE*0.05);
    thput_per_UE_5 = UE_thput_sort(index_5);
end
thput_per_cell_5 = thput_per_UE_5*n_UE_served_per_BS;

if SYS_config.photo
    [range,output] = NR_plot_CDF(UE_thput);
    figure;
    subplot(1,2,1);
    plot(range,output);
    xlabel('Throughput (bps)');
    ylabel('F(x)');
    title('asynchronization UE to UE Throughput CDF');
end

% 吞吐量和SINR关系
if SYS_config.photo
    UE_SINR = sort(UE_SINR);
    UE_thput = sort(UE_thput);
    subplot(1,2,2);
    plot(UE_SINR,UE_thput);
    xlabel('SINR (dB)');
    ylabel('Throughput (bps)');
    title('asynchronization UE to UE Throughput VS SINR');
end

% 上行
ul_UE_SINR = ul_all_SINR(activate_UE_id);
for i=1:length(ul_UE_SINR)
    if ul_UE_SINR(i)<-12
        ul_UE_thput(i) = 0;
    elseif ul_UE_SINR(i)>22
        ul_UE_thput(i) = 0.4 * bandwidth/1e6*log2(1+10^(0.1*22));
    else
        ul_UE_thput(i) = 0.4 * bandwidth/1e6*log2(1+10^(0.1*ul_UE_SINR(i)));
    end
end

ul_thput_per_UE = sum(ul_UE_thput)/n_UE;
ul_thput_per_cell = ul_thput_per_UE*n_UE_served_per_BS;

if n_snapshot~=1
    ul_UE_thput_divi = zeros(n_snapshot,n_UE/n_snapshot);
    for i = 1:n_snapshot
        index_UE = [];
        for a_ = 1:n_UE_served_per_BS
            if isempty(index_UE)
                index_UE = (i*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE;
            else
                index_UE = [index_UE,(i*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE ];
            end
        end
        ul_UE_thput_divi(i,:) = ul_UE_thput(index_UE);
    end
    ul_UE_thput_divi = sort(ul_UE_thput_divi,2);
    ul_UE_thput_ave = sum(ul_UE_thput_divi)/n_snapshot; %得到平均吞吐量
    
    index_5 = ceil(n_UE/n_snapshot*0.05);
    ul_thput_per_UE_5 = ul_UE_thput_ave(index_5);
else
    ul_UE_thput_sort = sort(ul_UE_thput);
    index_5 = ceil(n_UE*0.05); %得到5%-tile吞吐量
    ul_thput_per_UE_5 = ul_UE_thput_sort(index_5);
end
ul_thput_per_cell_5 = ul_thput_per_UE_5*n_UE_served_per_BS;

if SYS_config.photo
    [range,output] = NR_plot_CDF(ul_UE_thput);
    subplot(1,2,2);
    plot(range,output);
    xlabel('Throughput (bps)');
    ylabel('F(x)');
    title('UL Throughput CDF');
end
if SYS_config.photo
    ul_UE_SINR = sort(ul_UE_SINR);
    ul_UE_thput = sort(ul_UE_thput);
    subplot(1,2,2);
    plot(ul_UE_SINR,ul_UE_thput);
    xlabel('SINR (dB)');
    ylabel('Throughput (bps)');
    title('UL Throughput VS SINR');
end

fprintf('UE to UE average per cell: %f UE to UE edge per cell: %f \n',thput_per_cell,thput_per_cell_5);
fprintf('BS to BS average per cell: %f BS to BS edge per UE: %f \n',ul_thput_per_cell,ul_thput_per_cell_5);
end