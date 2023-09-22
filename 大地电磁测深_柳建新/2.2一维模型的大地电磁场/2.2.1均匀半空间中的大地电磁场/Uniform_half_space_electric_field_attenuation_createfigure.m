function Uniform_half_space_electric_field_attenuation_createfigure(X1, Y1, X2, X3, X4)
%CREATEFIGURE(X1, Y1, X2, X3, X4)
%  X1:  x ���ݵ�����
%  Y1:  y ���ݵ�����
%  X2:  x ���ݵ�����
%  X3:  x ���ݵ�����
%  X4:  x ���ݵ�����

%  �� MATLAB �� 13-Sep-2023 17:38:46 �Զ�����

% ���� figure
figure1 = figure;

% ���� subplot
subplot1 = subplot(1,2,1,'Parent',figure1);
hold(subplot1,'on');

% ���� plot
plot(X1,Y1,'Parent',subplot1,'DisplayName','real','LineWidth',1.5);

% ���� plot
plot(X2,Y1,'Parent',subplot1,'DisplayName','imag','LineWidth',1.5,...
    'LineStyle','--');

% ���� ylabel
ylabel('Depth/km','FontSize',14,'FontName','Times New Roman');

% ���� xlabel
xlabel('Ex/Ex0','FontSize',14,'FontName','Times New Roman');

box(subplot1,'on');
% ������������������
set(subplot1,'FontName','Times New Roman');
% ���� legend
legend1 = legend(subplot1,'show');
set(legend1,'Location','southeast');

% ���� subplot
subplot2 = subplot(1,2,2,'Parent',figure1);
hold(subplot2,'on');

% ���� plot
plot(X3,Y1,'Parent',subplot2,'DisplayName','real','LineWidth',1.5);

% ���� plot
plot(X4,Y1,'Parent',subplot2,'DisplayName','imag','LineWidth',1.5,...
    'LineStyle','--');

% ���� ylabel
ylabel('Depth/km','FontSize',14,'FontName','Times New Roman');

% ���� xlabel
xlabel('Ex/Ex0','FontSize',14,'FontName','Times New Roman');

box(subplot2,'on');
% ������������������
set(subplot2,'FontName','Times New Roman');
% ���� legend
legend2 = legend(subplot2,'show');
set(legend2,'Location','southeast');

