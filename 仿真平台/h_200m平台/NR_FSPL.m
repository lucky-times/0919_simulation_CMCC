function PL = NR_FSPL(d,f)
    %   d:三维距离（m）
    %   f:载波频率（Hz）
    %   38.901 UMa LOS情况下的路损模型
    PL = 28 + 22*log10(d) + 20*log10(f/1e9);
%     PL = 32.45 + 20*log10(d) + 20*log10(f/1e9);
end

