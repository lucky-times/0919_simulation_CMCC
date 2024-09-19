function P = NR_randperm(N,K)
% 低版本的matlab不支持randperm(N,K)这种样式的输入
% 输入参数
% N：产生一个1到N的整数序列，并打乱
% K：从上面序列中取前K个
% 输出参数
% P：一个打乱顺序了的1到N的整数序列，取前K个

P_1 = randperm(N);
try
    P = P_1(1:K);
catch
    % 避免出现K大于N的情况
    error('K must <= N');
end

end

