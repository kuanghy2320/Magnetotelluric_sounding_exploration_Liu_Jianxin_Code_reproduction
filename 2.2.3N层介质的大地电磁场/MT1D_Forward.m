%%
warning off
clc
clear all

%%%%%%%%%%%%%%%%%%%%%%
%%����
% rho=[1 0.000000000000000001;1 0.0025;1 0.025;1 0.1;1 0.25;1 1;1 3;1 10;1 25;1 200;1 10000000000000000000000];
% h=[20000];
% 
% figure(1)
% for i=1:size(rho,1)
%     [b,T,y,pc,ph]=MT1D_Forward1(rho(i,:),h);
%     loglog(T,pc./rho(1));
%     hold on
% end
% hold off
% figure(2)
% for i=1:size(rho,1)
%     [b,T,y,pc,ph]=MT1D_Forward1(rho(i,:),h);
%     semilogx(T,ph);
%     hold on
% end
% hold off
%%%%%%%%%%%%%%%%%%%%%%

%%
%%����
rho=[100 500 3000];
h=[100 500;100 2000;100 6000;100 10000;100 20000;100 50000]
figure(1)
for i=1:size(h,1)
[b,T,y,pc,ph]=MT1D_Forward1(rho,h(i,:));
loglog(T,pc./rho(1));
hold on
end
hold off
figure(2)
for i=1:size(h,1)
[b,T,y,pc,ph]=MT1D_Forward1(rho,h(i,:));
semilogx(T,ph);
hold on
end
hold off
%%
function [b,T,y,pc,ph]=MT1D_Forward1(rho,h)
mu=(4e-7)*pi;
T=logspace(-5,7,100);
omega=2*pi./T;
i=sqrt(-1);
k=zeros(size(rho,2),size(T,2));
for N=1:size(rho,2)
    k(N,:)=sqrt(-i*2*pi*mu./(T.*rho(N)));%ÿһ�㸴����
end

m=size(rho,2);
%%�迹���㹫ʽ
y=-i*2*pi*mu./(T.*k(m,:));
for nn=m-1:-1:1
    A=-i*2*pi*mu./(T.*k(nn,:));
    B=exp(-2*k(nn,:).*h(nn));
    y=A.*(A.*(1-B)+y.*(1+B))./(A.*(1+B)+y.*(1-B));
end
pc=abs(y).^2./(mu*2*pi./T);
ph=-atan(imag(y)./real(y)).*180/pi;
b=[pc,ph];
end

% function createfigure1(X1, YMatrix1)
% %CREATEFIGURE(X1, YMatrix1)
% %  X1:  x ���ݵ�����
% %  YMATRIX1:  y ���ݵľ���
% 
% %  �� MATLAB �� 13-Oct-2023 14:51:20 �Զ�����
% 
% % ���� figure
% figure1 = figure;
% 
% % ���� axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% 
% % ʹ�� loglog �ľ������봴������
% loglog(X1,YMatrix1,'Parent',axes1);
% 
% % ���� ylabel
% ylabel('Pa/p1');
% 
% % ���� xlabel
% xlabel('T/s');
% 
% % ���� title
% title('��������ӵ�������Ӧ����');
% 
% % ȡ�������е�ע���Ա����������� X ��Χ
% % xlim(axes1,[100 10000000000]);
% box(axes1,'on');
% % ������������������
% set(axes1,'FontName','����','FontSize',12,'XMinorTick','on','XScale','log',...
%     'YMinorTick','on','YScale','log');
% end
