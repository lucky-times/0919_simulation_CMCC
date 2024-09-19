function networkShadowFadingMap = NR_generate_shadowfading(networkPathlossMap,eNodeB_sites,SYS_config)
%用于产生阴影衰落
%输入参数：
%networkPathlossMap：路损图
%eNodeBs：eNodeB实体
%SYS_config：参数配置
%
%输出参数：
%networkShadowFadingMap：阴影衰落图谱
if SYS_config.debug_level>=1
    fprintf('Generating shadow fading\n');
end

networkShadowFadingMap = channel_gain_wrappers.shadowFadingMapClaussen(SYS_config,networkPathlossMap,eNodeB_sites);

if SYS_config.debug_level>=2
    networkShadowFadingMap.print;
end

end