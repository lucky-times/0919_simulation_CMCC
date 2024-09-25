clear all;
close all;
clc
%% 100m
r = [100, 300, 500, 800, 1200];
mean_overlap_rate = zeros(1,5);
for j = 1:length(r)
    str1 = sprintf("PCI_%dm_Measured.csv", r(j));
    str2 = sprintf("PCI_%dm_Simulated.csv", r(j));
    tmp1 = readmatrix(str1);
    tmp2 = readmatrix(str2);
    row = length(tmp1);
    overlap_rate = zeros(row, 1);
    for i = 1:row
        tmp = tmp1(i, :);
        num = sum(~isnan(tmp));
        if num <= 8
            % 找出交集
            C = intersect(tmp(~isnan(tmp)), tmp2(i, :));        
            % 计算交集的个数
            numIntersections = numel(C);
            overlap_rate(i) = numIntersections/num;
        else
            % 找出交集
            C = intersect(tmp(~isnan(tmp)), tmp2(i, :));        
            % 计算交集的个数
            numIntersections = numel(C);
            if numIntersections <= 8
                overlap_rate(i) = numIntersections/num;
            else
                overlap_rate(i) = 1;
            end
        end
    end
    mean_overlap_rate(j) = mean(overlap_rate);
    str3 = sprintf('PCI_%dm_overlap_rate.csv', r(j));
    writematrix(overlap_rate, str3);
    % 打开文件以追加内容
    fileID = fopen(str3, 'a');  % 'a' 表示追加模式
    % 检查文件是否成功打开
    if fileID == -1
        error('File could not be opened for appending.');
    end
    % 定义要追加的字符串
    strToAdd = sprintf('预测干扰小区PCI重合率平均值为：%.2f', mean_overlap_rate(j));
    % 将字符串写入文件的最后一行下面
    fprintf(fileID, strToAdd);
    % 关闭文件
    fclose(fileID);
end