clc
clear all

%%
mu=4*pi*1e-7;
DZ=1;%��Ԫ����
NZ=100000;%��Ԫ����
rho=ones(1,NZ);%��Ԫ�����ʣ����Ȱ�ռ�


%f=10.^[-2 -1 0 1 2 3];%Ƶ��
f=logspace(-2,3,100);
Depth=-DZ/2:-DZ:-DZ*NZ-DZ/2;%���
% depth=[0.5,1,2,4,,30000]
z=-Depth;
t=0;%ʱ��

[pc,ph,Ex]=MT1D_Forward_FDM1(rho,mu,f,NZ,DZ,Depth);%���޲�ּ���

%%
figure(1)
loglog(f,pc);

figure(2)
plot(f,ph)

figure(3)
%���ۼ���
p=rho(1);
Ex_Ex0=exp(-i.*(2*pi*f.*t-sqrt(2*pi*f.*mu/p/2).*z')-sqrt(2*pi*f.*mu/p/2).*z');

x=Depth(1:end-1);
y1r=real(Ex_Ex0(1:end-1,50));y1i=imag(Ex_Ex0(1:end-1,50));
y2r=real(Ex(:,50));y2i=imag(Ex(:,50));
plot(x,y1r,x,y1i);
hold on
plot(x,y2r,x,y2i);

%������
figure(4)
plot(x,(y2r-y1r)./y1r,x,(y2i-y1i)./y1i)

%%
function [pc,ph,Ex]=MT1D_Forward_FDM1(rho,mu,f,NZ,DZ,Depth)
    %rhoΪ�ڵ㴦�ĵ�����ֵ
    %NZΪ��Ԫ����
    %DZΪ��Ԫ����
    eps=8.8419e-012;
    
    for ff=1:size(f,2)
        omega=2*pi*f(ff);
        
        P=sparse(NZ,1);%ȫ��ϡ�����
        K=sparse(NZ,NZ);%ȫ��ϡ�����
        K(1,1:2)=[1,1];
        P(1)=2;
        for j=2:NZ-1
            K(j,j-1:j+1)=[1/DZ^2,(0.+1i)*omega*mu/rho(j)-2/DZ^2,1/DZ^2];
        end
        K(NZ,NZ)=1;
        Ex(:,ff)=K\P;%���Է�������⣺ֱ�ӷ�

        Ex_g=(Ex(1,ff)+Ex(2,ff))/2;
        Hy_g=(Ex(2,ff)-Ex(1,ff))/DZ/((0.+1i)*mu*omega);
        
        Z(ff)=Ex_g/Hy_g;
        pc(ff)=abs(Z(ff))^2/mu/omega;
        ph(ff)=-atan(imag( Z(ff))/real( Z(ff)))*180/pi;
    end
end