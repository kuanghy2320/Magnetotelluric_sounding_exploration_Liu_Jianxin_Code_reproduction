clc
clear all

DZ=20.0;%单元长度
NZ=100;%单元个数
rho=ones(1,4+(NZ-1)*3);%节点电阻率，均匀半空间
Depth=0:-DZ/3:-DZ*NZ;%深度
%f=10.^[0 1 2];%频率
f=logspace(-5,7,100)
t=0;
mu=(4e-7)*pi;

[pc,ph,Ex]=MT1D_Forward_FEM1(rho,f,NZ,DZ)
%%
figure(1)
loglog(f,pc);

figure(2)
plot(f,ph)

figure(3)
p=rho(1);
z=-Depth;
Ex_Ex0=exp(-i.*(2*pi*f.*t-sqrt(2*pi*f.*mu/p/2).*z')-sqrt(2*pi*f.*mu/p/2).*z');

plot(Depth,real(Ex_Ex0(:,3)),Depth,imag(Ex_Ex0(:,3)))
hold on
plot(Depth,real(Ex(:,3)),Depth,-imag(Ex(:,3)))

%%
function [pc,ph,Ex]=MT1D_Forward_FEM1(rho,f,NZ,DZ)
    %rho为节点处的电阻率值
    %NZ为单元个数
    %DZ为单元长度
    
    NE=NZ;%单元总数
    NP=4+(NE-1)*3;%节点总数
    L=ones(1,NE)*DZ;%单元长度
    
    %%存放单元节点编号
    for i=1:NE
        ME(1,i)=1+(i-1)*3;
        ME(2,i)=2+(i-1)*3;
        ME(3,i)=3+(i-1)*3;
        ME(4,i)=4+(i-1)*3;
    end

    for ff=1:size(f,2)
        K1=sparse(NP,NP);K2=sparse(NP,NP);K3=sparse(NP,NP);P=sparse(NP,1);
%%
%%%%%%%%%%%%% 形成K1e矩阵 %%%%%%%%%%%%%%%
        %遍历所有单元
        for h=1:NE
            l=L(h);%l表示单元长度
            K=[148 -189 54 -13;-189 432 -297 54;54 -297 432 -189;-13 54 -189 148];
            for j=1:4
                NJ=ME(j,h);
                for k=1:4
                    NK=ME(k,h);
                    K1(NJ,NK)=K1(NJ,NK)+K(j,k)/40/l;
                end
            end
        end
%%
%%%%%%%%%%%%% 形成K2e矩阵 %%%%%%%%%%%%%%%
        %遍历所有单元
        mu=(4e-7)*pi;
        w=2*pi*f(ff);
        m=sqrt(-1)*w*mu;
        for h=1:NE
            l=L(h);
            r1=1/rho(1+(h-1)*3);%节点电导率值
            r2=1/rho(2+(h-1)*3);%节点电导率值
            r3=1/rho(3+(h-1)*3);%节点电导率值
            r4=1/rho(4+(h-1)*3);%节点电导率值
            K(1,1)=356*r1+  216*r2-  81*r3+ 20*r4;
            K(2,1)=216*r1+  324*r2-162*r3+ 18*r4;
            K(2,2)=356*r1+  216*r2-    0*r3+ 81*r4;
            K(3,1)=-81*r1 -  162*r2+  81*r3+ 18*r4;
            K(3,2)=-162*r1+216*r2-     0*r3-162*r4;
            K(3,3)=   81*r1+    0*r2+2187*r3+324*r4;
            K(4,1)=20*r1+    18*r2+18*r3+20*r4;
            K(4,2)=18*r1+81*r2-162*r3-81*r4;
            K(4,3)=18*r1-162*r2+324*r3+216*r4;
            K(4,4)=20*r1-81*r2+216*r3+357*r4;
            
            K(2,1)=K(1,2);K(3,1)=K(1,3);K(4,1)=K(1,4);
                                K(3,2)=K(2,3);K(4,2)=K(2,4);
                                                    K(4,3)=K(3,4);
            for j=1:4
                NJ=ME(j,h);
                for k=1:4
                    NK=ME(k,h);
                    K1(NJ,NK)=K1(NJ,NK)+K(j,k)*m*l/6720;
                end
            end
        end
%%
%%%%%%%%%%%%% 形成K3e矩阵 %%%%%%%%%%%%%%%
        mu=(4e-7)*pi;
        w=2*pi*f(ff);
        m=sqrt(-1)*w*mu;
        a=sqrt(-m/rho(NP));
        K3(NP,NP)=a;
        
%%
%%%%%%%%%%%%% 组装总体刚度矩阵 %%%%%%%%%%%%%%%
        v=K1-K2+K3;%总体合成
        v(NP,:)=0;
        v(NP,NP)=1;
        
%%
%%%%%%%%%%%%% 代入边值条件 %%%%%%%%%%%%%%%
        DJ=1;%地表节点
        v(DJ,DJ)=v(DJ,DJ)*10^10;
        P(DJ)=v(DJ,DJ)*1;%上边界条件
        
%%
%%%%%%%%%%%%% 线性方程组直接解法 %%%%%%%%%%%%%%%
        Ex(:,ff)=v\P;%电场Ex
        
        DE(ff)=(-11*Ex(DJ,ff)+18*Ex(DJ+1,ff)-9*Ex(DJ+2,ff)+2*Ex(DJ+3,ff))/L(11)/2;
        Z(ff)=DE(ff)/(sqrt(-1)*w*mu);
        pc(ff)=abs(Z(ff))^2/mu/w;
        ph(ff)=-atan(imag( Z(ff))/real( Z(ff)))*180/pi;
    end
    
end