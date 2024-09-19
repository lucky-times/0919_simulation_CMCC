data2 = readtable('7个4.9G站点.xlsx');
% 经度、纬度
info = data2(:, [9, 10]);
% 小区标识、方位角
info2 = data2(:, [23, 31]);
eNodeB_height = data2(:, 32);
info = table2array(info);
pos = [];
for i = 1:length(info)
    if(ismember(info(i, 1), pos))
        continue;
    else
        pos = [pos;info(i, :)];
    end
end


% %% 将csv文件转化为excel，便于观察数据信息
% % 指定CSV文件的路径
% csvFilePath = 'D:\\OneDrive\\OneDrive - bupt.edu.cn\\桌面\\中国移动研究院\\外场测试数据\7.4.2多层窄波束配置--蜂窝三扇区组网\\7.4.2多层窄波束配置--蜂窝三扇区组网-1km站间距-飞行高度200米-加扰50%.csv';
% 
% % 读取CSV文件
% data = readtable(csvFilePath);
% 
% % 指定输出的Excel文件路径
% excelFilePath = 'D:\OneDrive\OneDrive - bupt.edu.cn\桌面\中国移动研究院\外场测试数据\7.4.2多层窄波束配置--蜂窝三扇区组网\7.4.2多层窄波束配置--蜂窝三扇区组网-1km站间距-飞行高度200米-加扰50%.xlsx';
% 
% % 写入Excel文件
% writetable(data, excelFilePath);

%% 测试
k1 = -174;
B = 20e6;
k1 + 10*log10(B)

