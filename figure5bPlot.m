clear;
clc
close all
% clip level 2
load('figure5a.mat')
plotimage(Vx2)
xlabel('x/dx')
ylabel('z/dz')
title('')
aa=Vx2;


load('figure5b.mat')
plotimage(Vx2)
xlabel('x/dx')
ylabel('z/dz')
title('')
bb=Vx2;

load('figure5c.mat')
plotimage(Vx2)
xlabel('x/dx')
ylabel('z/dz')
title('')
cc=Vx2;


figure;plot(aa(:,400),'b','linewidth',2)
hold on;plot(bb(:,400),'r','linewidth',2);
hold on;plot(cc(:,400),'k','linewidth',2);



legend('Previous Non-balanced scheme:A','Improved Non-balanced scheme:B','Tan & Huang:C')
axis([0 565 -1.7*10^-7 2.5*10^-7])
grid on
xlabel('z/dz');ylabel('Amp')