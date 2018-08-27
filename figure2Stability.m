clear;
clc
close all
global M dt h v;

v=1500;
h=15;
for loop=1:3
    if loop==1
        M=3;
        k=linspace(1/50,0.35*pi/h,M);
        iii=1;
        for dt=0.00005:0.00005:0.01;
            
            rr(iii)=v*dt/h;
            x0=0.001*ones(1,M+1);
            options = optimset('TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',8000,'MaxIter',200);
            
            lb=-5*ones(M+1,1);
            ub=5*ones(M+1,1);
            [c,fval,out,iteration]= fmincon(@myfun,x0,[],[],[],[],lb,ub,[],options) ;   % Invoke optimizer
            
            temp=0;
            for jjj=1:M
                temp=(-1)^jjj*c(jjj)*4+temp;
            end
            temp=temp+8*c(M+1);
            ss(iii)=sqrt(-2/temp);
            iii=iii+1;
        end
        figure; plot(rr,(rr),'r','linewidth',2)
        hold on;plot(rr,(ss),'k','linewidth',2);
    elseif loop==2
        M=5;
        k=linspace(1/50,0.55*pi/h,M);
        
        iii=1;
        for dt=0.00005:0.00005:0.01;
            
            rr(iii)=v*dt/h;
            x0=0.001*ones(1,M+1);
            options = optimset('TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',8000,'MaxIter',200);
            
            lb=-5*ones(M+1,1);
            ub=5*ones(M+1,1);
            [c,fval,out,iteration]= fmincon(@myfun,x0,[],[],[],[],lb,ub,[],options) ;   % Invoke optimizer
            
            temp=0;
            for jjj=1:M
                temp=(-1)^jjj*c(jjj)*4+temp;
            end
            temp=temp+8*c(M+1);
            ss(iii)=sqrt(-2/temp);
            iii=iii+1;
        end
        
        hold on;plot(rr,(ss),'b','linewidth',2);
    else
        M=7;
        k=linspace(1/50,0.78*pi/h,M);
        iii=1;
        for dt=0.00005:0.00005:0.01;
            
            rr(iii)=v*dt/h;
            x0=0.001*ones(1,M+1);
            options = optimset('TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',8000,'MaxIter',200);
            
            lb=-5*ones(M+1,1);
            ub=5*ones(M+1,1);
            [c,fval,out,iteration]= fmincon(@myfun,x0,[],[],[],[],lb,ub,[],options) ;   % Invoke optimizer
            
            temp=0;
            for jjj=1:M
                temp=(-1)^jjj*c(jjj)*4+temp;
            end
            temp=temp+8*c(M+1);
            ss(iii)=sqrt(-2/temp);
            iii=iii+1;
        end
        hold on;plot(rr,(ss),'m','linewidth',2);
    end
    
    
end
legend('r','M=3','M=5','M=7')
xlabel('r')
ylabel('Stability')
grid on
