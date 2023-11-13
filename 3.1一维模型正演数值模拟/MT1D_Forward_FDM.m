clc
clear all

%%
mu=4*pi*1e-7;
DZ=1;%单元长度
NZ=100000;%单元个数
rho=ones(1,NZ);%单元电阻率，均匀半空间


%f=10.^[-2 -1 0 1 2 3];%频率
f=logspace(-2,3,100);
Depth=-DZ/2:-DZ:-DZ*NZ-DZ/2;%深度
% depth=[0.5,1,2,4,,30000]
z=-Depth;
t=0;%时间

[pc,ph,Ex]=MT1D_Forward_FDM1(rho,mu,f,NZ,DZ,Depth);%有限差分计算

%%
figure(1)
loglog(f,pc);

figure(2)
plot(f,ph)

figure(3)
%理论计算
p=rho(1);
Ex_Ex0=exp(-i.*(2*pi*f.*t-sqrt(2*pi*f.*mu/p/2).*z')-sqrt(2*pi*f.*mu/p/2).*z');

x=Depth(1:end-1);
y1r=real(Ex_Ex0(1:end-1,50));y1i=imag(Ex_Ex0(1:end-1,50));
y2r=real(Ex(:,50));y2i=imag(Ex(:,50));
plot(x,y1r,x,y1i);
hold on
plot(x,y2r,x,y2i);

%误差检验
figure(4)
plot(x,(y2r-y1r)./y1r,x,(y2i-y1i)./y1i)

%%
function [pc,ph,Ex]=MT1D_Forward_FDM1(rho,mu,f,NZ,DZ,Depth)
    %rho为节点处的电阻率值
    %NZ为单元个数
    %DZ为单元长度
    eps=8.8419e-012;
    
    for ff=1:size(f,2)
        omega=2*pi*f(ff);
        
        P=sparse(NZ,1);%全零稀疏矩阵
        K=sparse(NZ,NZ);%全零稀疏矩阵
        K(1,1:2)=[1,1];
        P(1)=2;
        for j=2:NZ-1
            K(j,j-1:j+1)=[1/DZ^2,(0.+1i)*omega*mu/rho(j)-2/DZ^2,1/DZ^2];
        end
        K(NZ,NZ)=1;
        Ex(:,ff)=K\P;%线性方程组求解：直接法

        Ex_g=(Ex(1,ff)+Ex(2,ff))/2;
        Hy_g=(Ex(2,ff)-Ex(1,ff))/DZ/((0.+1i)*mu*omega);
        
        Z(ff)=Ex_g/Hy_g;
        pc(ff)=abs(Z(ff))^2/mu/omega;
        ph(ff)=-atan(imag( Z(ff))/real( Z(ff)))*180/pi;
    end
end