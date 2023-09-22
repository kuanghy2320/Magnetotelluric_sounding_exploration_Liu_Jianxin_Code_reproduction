function Uniform_half_space_electric_field_attenuation_createfigure(X1, Y1, X2, X3, X4)
%CREATEFIGURE(X1, Y1, X2, X3, X4)
%  X1:  x 数据的向量
%  Y1:  y 数据的向量
%  X2:  x 数据的向量
%  X3:  x 数据的向量
%  X4:  x 数据的向量

%  由 MATLAB 于 13-Sep-2023 17:38:46 自动生成

% 创建 figure
figure1 = figure;

% 创建 subplot
subplot1 = subplot(1,2,1,'Parent',figure1);
hold(subplot1,'on');

% 创建 plot
plot(X1,Y1,'Parent',subplot1,'DisplayName','real','LineWidth',1.5);

% 创建 plot
plot(X2,Y1,'Parent',subplot1,'DisplayName','imag','LineWidth',1.5,...
    'LineStyle','--');

% 创建 ylabel
ylabel('Depth/km','FontSize',14,'FontName','Times New Roman');

% 创建 xlabel
xlabel('Ex/Ex0','FontSize',14,'FontName','Times New Roman');

box(subplot1,'on');
% 设置其余坐标区属性
set(subplot1,'FontName','Times New Roman');
% 创建 legend
legend1 = legend(subplot1,'show');
set(legend1,'Location','southeast');

% 创建 subplot
subplot2 = subplot(1,2,2,'Parent',figure1);
hold(subplot2,'on');

% 创建 plot
plot(X3,Y1,'Parent',subplot2,'DisplayName','real','LineWidth',1.5);

% 创建 plot
plot(X4,Y1,'Parent',subplot2,'DisplayName','imag','LineWidth',1.5,...
    'LineStyle','--');

% 创建 ylabel
ylabel('Depth/km','FontSize',14,'FontName','Times New Roman');

% 创建 xlabel
xlabel('Ex/Ex0','FontSize',14,'FontName','Times New Roman');

box(subplot2,'on');
% 设置其余坐标区属性
set(subplot2,'FontName','Times New Roman');
% 创建 legend
legend2 = legend(subplot2,'show');
set(legend2,'Location','southeast');

