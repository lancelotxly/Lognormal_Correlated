clear all;clc;close all;

%%% PDF and intergral region transimit%%
% p=0.99
% g_1=-6:0.1:6;
% g_2=-6:0.1:6;
% [G_1,G_2]=meshgrid(g_1,g_2);
% u_G=2;
% sigma_G=1.2;
% part1=1/sqrt((2*pi)^2*((sigma_G^2)^2-(p*sigma_G^2)^2));
% part2=((G_1-u_G).^2+(G_2-u_G).^2-2*p.*(G_1-u_G).*(G_2-u_G))./(2*sigma_G^2*(p^2-1));
% f=part1.*exp(part2);
% contour(G_1,G_2,f,'r');hold on;
% h1=ezplot('exp(x1)+exp(x2)-sqrt(2)');hold on;
% set(h1,'Color','r');
% 
% p=0;
% g_1=-6:0.1:6;
% g_2=-6:0.1:6;
% [G_1,G_2]=meshgrid(g_1,g_2);
% u_G=2;
% sigma_G=1.2;
% part1=1/sqrt((2*pi)^2*((sigma_G^2)^2-(p*sigma_G^2)^2));
% part2=((G_1-u_G).^2+(G_2-u_G).^2-2*p.*(G_1-u_G).*(G_2-u_G))./(2*sigma_G^2*(p^2-1));
% f=part1.*exp(part2);
% contour(G_1,G_2,f,'b');hold on;
%  
% p=0.99;
% a=(1-sqrt(1-p^2))/p;
% h1=ezplot('exp(0.867608727478122*x1+x2)+exp(x1+0.867608727478122*x2)-sqrt(2)');hold on;
% set(h1,'Color','b');
% 
% x1=-6:0.1:6;
% xm_1=1/(a+1)*log(sqrt(1/2));
% xm_2=xm_1;
% k=(a*exp(a*xm_1+xm_2)+exp(xm_1+a*xm_2))/(exp(a*xm_1+xm_2)+a*exp(xm_1+a*xm_2));
% x2=-k.*(x1-xm_1)+xm_2;
% plot(x1,x2,'k--');hold on;
% 
% x1=[-6,xm_1];
% x2=[xm_2,xm_2];
% plot(x1,x2,'k--');
% x1=[xm_1,xm_1];
% x2=[-6,xm_2];
% plot(x1,x2,'k--');hold on;

%%%%%% only integral region transimit%%
% p=0.01 transimited region
p1=0.01;
a=(1-sqrt(1-p1^2))/p1;
h_p1=ezplot('exp(0.005000125006249*x1+x2)+exp(x1+0.005000125006249*x2)-sqrt(2)');hold on;
set(h_p1,'Color','b');
x1=-6:0.1:6;
xm_1=1/(a+1)*log(sqrt(1/2));
xm_2=xm_1;
k=(a*exp(a*xm_1+xm_2)+exp(xm_1+a*xm_2))/(exp(a*xm_1+xm_2)+a*exp(xm_1+a*xm_2));
x2=-k.*(x1-xm_1)+xm_2;
b_1=plot(x1,x2,'b--');hold on;
x1=[-6,xm_1];
x2=[xm_2,xm_2];
plot(x1,x2,'b--');
x1=[xm_1,xm_1];
x2=[-6,xm_2];
plot(x1,x2,'b--');hold on;

% p=0.6 transimited region
p2=0.6;
a=(1-sqrt(1-p2^2))/p2;
h_p2=ezplot('exp(00.333333333333333*x1+x2)+exp(x1+0.333333333333333*x2)-sqrt(2)');hold on;
set(h_p2,'Color','g');
x1=-6:0.1:6;
xm_1=1/(a+1)*log(sqrt(1/2));
xm_2=xm_1;
k=(a*exp(a*xm_1+xm_2)+exp(xm_1+a*xm_2))/(exp(a*xm_1+xm_2)+a*exp(xm_1+a*xm_2));
x2=-k.*(x1-xm_1)+xm_2;
b_2=plot(x1,x2,'g--');hold on;
x1=[-6,xm_1];
x2=[xm_2,xm_2];
plot(x1,x2,'g--');
x1=[xm_1,xm_1];
x2=[-6,xm_2];
plot(x1,x2,'g--');hold on;

% p=0.99 transimited region
p3=0.99;
a=(1-sqrt(1-p3^2))/p3;
h_p3=ezplot('exp(0.867608727478122*x1+x2)+exp(x1+0.867608727478122*x2)-sqrt(2)');hold on;
set(h_p3,'Color','r');
x1=-6:0.1:6;
xm_1=1/(a+1)*log(sqrt(1/2));
xm_2=xm_1;
k=(a*exp(a*xm_1+xm_2)+exp(xm_1+a*xm_2))/(exp(a*xm_1+xm_2)+a*exp(xm_1+a*xm_2));
x2=-k.*(x1-xm_1)+xm_2;
b_3=plot(x1,x2,'r--');hold on;
x1=[-6,xm_1];
x2=[xm_2,xm_2];
plot(x1,x2,'r--');
x1=[xm_1,xm_1];
x2=[-6,xm_2];
plot(x1,x2,'r--');hold on;

grid on;
grid minor;
legend([h_p1,b_1,h_p2,b_2,h_p3,b_3],'\Phi(\gamma_{th}),\rho=0.01','\Phi^{+}(\gamma_{th}) or \Phi^{-}(\gamma_{th}),\rho=0.01',...
    '\Phi(\gamma_{th}),\rho=0.6','\Phi^{+}(\gamma_{th}) or \Phi^{-}(\gamma_{th}),\rho=0.6',...
    '\Phi(\gamma_{th}),\rho=0.99','\Phi^{+}(\gamma_{th}) or \Phi^{-}(\gamma_{th}),\rho=0.99',...
 'Location','northeast');

xlabel({'$X_{1}$'},'Interpreter','latex');
ylabel({'$X_{2}$'},'Interpreter','latex');
title('');