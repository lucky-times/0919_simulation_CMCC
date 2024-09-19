function [ pos ] = create_eNodeBs(ISD,n_rings)
        %创建基站位置
        
        [tmp_gridx,tmp_gridy] = meshgrid(-n_rings:n_rings,...
            (-n_rings:n_rings)*sin(pi/3));
        % 产生点，类似于下面这种形式
        % * * * * *
        % * * * * *
        % * * * * *
        % * * * * *
        % * * * * *
        if mod(n_rings,2) == 0
            tmp_shift_idx = 2:2:2*n_rings+1; %转换偶数列
            % 如下
            % * * * * *
            %  * * * * *
            % * * * * *
            %  * * * * *
            % * * * * *
        else
            tmp_shift_idx = 1:2:2*n_rings+1; %转换奇数列
            % 如下
            %  * * * * *
            % * * * * *
            %  * * * * *
            % * * * * *
            %  * * * * *
        end
        
        tmp_gridx(tmp_shift_idx,:) = tmp_gridx(tmp_shift_idx,:) + 0.5;%转换
        
        rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)]; %rot是一个函数，w_是其输入参数，函数的作用是产生一个矩阵 顺时针旋转w_
        for i_ = 1:7
            tmp_hex(i_,:) = ((n_rings+0.5)*rot(pi/3)^(i_-1)*[1;0]).'; % 乘[1;0]表示只取第一列
        end
        
        tmp_valid_positions = inpolygon(tmp_gridx,tmp_gridy,tmp_hex(:,1),tmp_hex(:,2));% tmp_hex是正六边形顶点坐标，inpolygon判断点是否在六边形内
        tmp_x = tmp_gridx(tmp_valid_positions);
        tmp_y = tmp_gridy(tmp_valid_positions);
        
        for b_ = 1:length(tmp_x)
            pos(:, b_) = [tmp_x(b_)*ISD tmp_y(b_)*ISD];
        end