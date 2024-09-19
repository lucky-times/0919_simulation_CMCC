function defaultAxes
%UNTITLED 设置坐标轴属性
%   
  ax = gca; hold on; box on
  ax.XGrid = 'on';
  ax.YGrid = 'on';
  ax.XMinorTick = 'on';
  ax.YMinorTick = 'on';
  ax.LineWidth = 1.2;
  ax.GridLineStyle = ':';
  ax.FontName = 'Cambria';
  ax.FontSize = 12;
  ax.GridAlpha = .5;
end