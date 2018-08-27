function F = myfun(x)

global M dt h v;
k=linspace(1/10000,0.8*pi/h,100);
r=v*dt/h;
F=zeros(1,100);

for jj=1:5
    
    xita=(jj-1)*pi/16;
    temp=0;
    for m=1:M
        temp=temp+2*x(m)*(cos(m*k*sin(xita)*h)-cos((m-1)*k*sin(xita)*h)+cos(m*k*cos(xita)*h)-cos((m-1)*k*cos(xita)*h));
    end
    temp=temp+x(M+1)*( 4*cos(k*sin(xita)*h).*(cos(k*h*cos(xita))-1) +4*cos(k*cos(xita)*h).*(cos(k*h*sin(xita))-1) );
    
    bb=2/r^2*(cos(k*v*dt)-1);
    F(jj,:)=temp./bb;
end

F=F-1;

F=F(1,:)*F(1,:)'+F(2,:)*F(2,:)'+F(3,:)*F(3,:)'+F(4,:)*F(4,:)'+F(5,:)*F(5,:)';