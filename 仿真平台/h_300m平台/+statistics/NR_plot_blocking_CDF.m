function NR_plot_blocking_CDF( SYS_config,interference_level )
% 绘制带内阻塞 cdf曲线
% 输入参数：
% SYS_config：LTE参数系统
% interference_level： 阻塞信号

if SYS_config.photo
    if SYS_config.isDouble
        % 绘制阻塞信号CDF
        interference_level_1 = interference_level;
        interference_level_1 = 10.*log10(interference_level_1.*1000);
        [range,output] = NR_plot_CDF(interference_level_1);
        
        % 为了能取到99.99%的点，对曲线进行平滑
        ft = fittype( 'smoothingspline' );
        opts = fitoptions( ft );
        opts.SmoothingParam = 1;
        [fitresult, ~] = fit( range', output', ft, opts );
        range = linspace(min(range),max(range),length(range)*100);
        output = feval(fitresult,range);
        
        figure;
        plot(range,output,'b');
        axis([-inf,inf,0,1]);
        grid on;
        xlabel('Blocking（dBm）');
        ylabel('F(x)');
        title('in-band blocking signal CDF');
        % 保存结果，为了将不同功控敏感度的结果画在同一个fig上
        file_name1 = ['.\tmp\','Blocking_range_',num2str((SYS_config.PCe_param)^2)];% 保存在tmp中，如果需要画不同功控错误情况下的图形，可以写个脚本调用
        file_name2 = ['.\tmp\','Blocking_output_',num2str((SYS_config.PCe_param)^2)];
        save(file_name1,'range');
        save(file_name2,'output');
    end
end


end

