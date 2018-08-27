clear;
clc
close all

load('figure5a.mat')
% plotimage(Vx2)
% xlabel('x/dx')
% ylabel('z/dz')
% title('')
aa=Vx2;


load('figure5b.mat')
% plotimage(Vx2)
% xlabel('x/dx')
% ylabel('z/dz')
% title('')
bb=Vx2;

load('figure5c.mat')
% plotimage(Vx2)
% xlabel('x/dx')
% ylabel('z/dz')
% title('')
cc=Vx2;


figure;plot(aa(:,400)-bb(:,400)+0.6*10^-7,'m-..','linewidth',2)
hold on;plot(bb(:,400)-cc(:,400)+0*10^-7,'g-..','linewidth',2)
hold on;plot(aa(:,400)-cc(:,400)-0.6*10^-7,'c-..','linewidth',2)


legend('A-B','B-C','A-C')
grid on
xlabel('z/dz')
ylabel('Difference')
axis([0 565 -1*10^-7 1.2*10^-7])