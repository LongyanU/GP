clear;
clc;
close all
global v dt h M
M=7;
% v=2800;
% dt=0.002;
% h=10;

dt=0.001;
v=2800;
h=20;

x0=0.001*ones(1,M+1);
options = optimset('TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',8000,'MaxIter',200);

lb=-5*ones(M+1,1);
ub=5*ones(M+1,1);
[x,fval,out,iteration]= fmincon(@myfun,x0,[],[],[],[],lb,ub,[],options)    % Invoke optimizer


r=v*dt/h;
k=linspace(1/1000000,pi/h,100);
F=zeros(5,400);

for kk=1:5
    
    xita=(kk-1)*pi/16;
    temp=0;
    for m=1:M
        temp=temp+x(m)*(cos(m*k*sin(xita)*h)-cos((m-1)*k*sin(xita)*h)+cos(m*k*cos(xita)*h)-cos((m-1)*k*cos(xita)*h));
    end
    temp=temp+1/2*x(M+1)*( 4*cos(k*sin(xita)*h).*(cos(k*h*cos(xita))-1) +4*cos(k*cos(xita)*h).*(cos(k*h*sin(xita))-1) );
    
    temp=1+temp*r^2;
    temp=acos(temp)./(k*v*dt);
    a1=(h/v*(1./temp-1));
    
    
    if kk==1
        figure;plot(a1,'k','LineWidth',2)
    elseif kk==2
        hold on;plot(a1,'m--','LineWidth',2);
    elseif kk==3
        hold on;plot(a1,'r:','LineWidth',2)
    elseif kk==4
        hold on; plot(a1,'b-.','LineWidth',2)
    else
        hold on;plot(a1,'c:.','LineWidth',2)
    end
    
end
legend('¦È=0', '¦È=¦Ð/16','¦È=2¦Ð/16','¦È=3¦Ð/16','¦È=4¦Ð/16')
xlabel('percentage of kh')
ylabel('\epsilon (\theta)')
grid on
axis([0 100 -1.5*10^-5  5*10^-5])