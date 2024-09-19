% 半径（自变量）
radii = [100, 300, 500, 800, 1200];

% 平均误差（因变量）
errors = [5.90537750549023, 6.02917497867582, 9.97924931062924, 2.93752556887793, 5.85669959810204];

% 使用 polyfit 进行线性回归拟合
[p, ~, mu] = polyfit(radii, errors, 2); % 1 表示线性模型（一次多项式）

% 使用拟合的模型来计算预测值
radii_linspace = linspace(min(radii), max(radii), 100); % 创建一个半径的线性空间用于绘图
errors_fit = polyval(p, radii_linspace, [], mu);

% 绘制原始数据和拟合曲线
figure;
plot(radii, errors, 'o', radii_linspace, errors_fit, '-');
xlabel('Radius (m)');
ylabel('Error');
title('Error vs. Radius Fit');
legend('Original Data', 'Fitted Curve');
grid on;