%%��������
f=[1 10] ;%��λHz
w=2*pi.*f;
u=4*pi*1e-7;%����еĴŵ���
p=10;%��λ��*m��������Ϊ10��*m�ľ��Ȱ�ռ�
t=0;%ʱ��
z=0:100:10000;%zΪ���
%%����
Ex_Ex0=exp(-i.*(w.*t-sqrt(w.*u/p/2).*z')-sqrt(w.*u/p/2).*z');

X1=real(Ex_Ex0(:,1));
Y1=-z/1000;
X2=imag(Ex_Ex0(:,1));
X3=real(Ex_Ex0(:,2));
X4=imag(Ex_Ex0(:,2));
Uniform_half_space_electric_field_attenuation_createfigure(X1, Y1, X2, X3, X4) 

% subplot(121)
% plot(real(Ex_Ex0(:,1)),-z/1000,"linewidth",1.5)
% xlabel("Ex/Ex0",'FontName','Times New Roman','FontSize',14)
% ylabel("Depth/km",'FontName','Times New Roman','FontSize',14)
% hold on
% plot(imag(Ex_Ex0(:,1)),-z/1000,'--',"linewidth",1.5)
% legend({"ʵ��ֵ","�鲿ֵ"},'Location','southeast')
% legend('boxoff')
% 
% subplot(122)
% plot(real(Ex_Ex0(:,2)),-z/1000,"linewidth",1.5)
% xlabel("Ex/Ex0",'FontName','Times New Roman','FontSize',14)
% ylabel("Depth/km",'FontName','Times New Roman','FontSize',14)
% hold on
% plot(imag(Ex_Ex0(:,2)),-z/1000,'--',"linewidth",1.5)
% legend({"ʵ��ֵ","�鲿ֵ"},'Location','southeast')
% legend('boxoff')