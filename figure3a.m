%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% authors Liang wenquan, Wang Yanfei,Yang changchun %%%%%%%%%%%%%
% Elapsed time is 20.682348 seconds.

clear
clc %%%%%%%
close all
nt=350;    % number of time steps
eps=.6;     % stability
isnap=60;    % snapshot sampling

nx=400;
nz=400;

v=ones(nz,nx)*2800;
v(1:nz/2,:)=2800;

p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=10;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.002; % calculate time step from stability criterion
tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0=160;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
zs=nz/2;
xs=nx/2;

h=dx;
r=v*dt/h;

for loop=1:1
    tic
    
    p=zeros([nz nx]); Vx=p; Vz=p;
    
    coeff=[ 1.3117, -0.154858, 0.0463066, -0.0165547, 0.00571777, -0.00159442, 0.000255275];
    for it=1:nt-2,
        
        d2px11=Vx-circshift(Vx,[0 1]);
        d2px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
        d2px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
        d2px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
        d2px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
        d2px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
        d2px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));
        
        %two point in Vz, six point in staggered grid FD
%         addPoint=(circshift(Vx,[-1 0])-circshift(Vx,[-1 1]));
%         addPoint=addPoint+(circshift(Vx,[1 0])-circshift(Vx,[1 1]));
        
        d2pz11=Vz-circshift(Vz,[1 0]);
        d2pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
        d2pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
        d2pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
        d2pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
        d2pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
        d2pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));
      
        d2px=coeff(1)*d2px11+coeff(2)*d2px12+coeff(3)*d2px13+coeff(4)*d2px14+coeff(5)*d2px15+coeff(6)*d2px16...
            +coeff(7)*d2px17;
        d2pz=coeff(1)*d2pz11+coeff(2)*d2pz12+coeff(3)*d2pz13+coeff(4)*d2pz14+coeff(5)*d2pz15+coeff(6)*d2pz16...
            +coeff(7)*d2pz17;
        
        p=p-dt*v.^2.*(d2px+d2pz)/h;
        p(zs,xs)= p(zs,xs)+src(it)*dt^2;
        
        
        d2px1=(circshift(p,[0 -1])-circshift(p,[0 0]));
        %         d2pz2=(circshift(p,[0 -2])-circshift(p,[0 1]));
        %         d2pz3=(circshift(p,[0 -3])-circshift(p,[0 2]));
        %         d2pz4=(circshift(p,[0 -4])-circshift(p,[0 3]));
        %         d2pz5=(circshift(p,[0 -5])-circshift(p,[0 4]));
        %         d2pz6=(circshift(p,[0 -6])-circshift(p,[0 5]));
        %         d2pz7=(circshift(p,[0 -7])-circshift(p,[0 6]));
        
        
        
        d2pz1=(circshift(p,[-1])-circshift(p,[0]));
        %         d2px2=(circshift(p,[-2])-circshift(p,[1]));
        %         d2px3=(circshift(p,[-3])-circshift(p,[2]));
        %         d2px4=(circshift(p,[-4])-circshift(p,[3]));
        %         d2px5=(circshift(p,[-5])-circshift(p,[4]));
        %         d2px6=(circshift(p,[-6])-circshift(p,[5]));
        %         d2px7=(circshift(p,[-7])-circshift(p,[6]));
        %            d2px=coeff(1)*d2px1+coeff(2)*d2px2+coeff(3)*d2px3+coeff(4)*d2px4+coeff(5)*d2px5+coeff(6)*d2px6...
        %             +coeff(7)*d2px7;
        %         d2pz=coeff(1)*d2pz1+coeff(2)*d2pz2+coeff(3)*d2pz3+coeff(4)*d2pz4+coeff(5)*d2pz5+coeff(6)*d2pz6...
        %             +coeff(7)*d2pz7;
        
        
        Vx=Vx-dt*d2px1/h;
        Vz=Vz-dt*d2pz1/h;
        
%         [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
        if rem(it,isnap)== 0,
            imagesc(x,z,p), axis equal
            colormap gray
            xlabel('x'),ylabel('z')
            title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
            drawnow
        end
        
        a(it)=p(nz/2-13,nx/2);
        b(it)=p(nz/2+13,nx/2);
        
        seis_record(it,:)=p(60,:);
        if it==500
            pp1=p;
        elseif it==1000
            pp2=p;
        elseif it==1500
            pp3=p;
        elseif it==2000
            pp4=p;
        elseif it==2500
            pp5=p;
        end
        
    end
    
end
toc


save('figure3a.mat')