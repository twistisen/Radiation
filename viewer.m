function viewer
close all
clc
A=load('spectrum.txt');
B=load('orbitsec0-0.txt');

pol1=0.5*(A(:,2)+A(:,3)+A(:,4)+A(:,5));
pol2=0.5*(A(:,6)+A(:,7)+A(:,8)+A(:,9));

E0=30e9;
m=511e3;
e=1/sqrt(137);
eta=1.0;
w0=1.55;
x=4*E0*w0/m^2;

figure(1)
plot(B(:,8),B(:,4))





f=figure(2)
set(f,'DefaultTextInterpreter','Latex')
set(f,'DefaultAxesTickLabelInterpreter','Latex')
#set(f,'DefaultlegendInterpreter','Latex')
plot(A(:,1)*1e-9,A(:,6))
hold all
plot(A(:,1)*1e-9,A(:,7))
plot(A(:,1)*1e-9,A(:,8))
% plot(A(:,1),A(:,5))

plot(A(:,1)*1e-9,A(:,2))
plot(A(:,1)*1e-9,A(:,3))
% plot(A(:,1),A(:,8),'-.')
% plot(A(:,1),A(:,9),'-.')

% plot(C(:,1)*1e-9,C(:,6),'black:','linewidth',1.5)
% hold all
% plot(C(:,1)*1e-9,C(:,7),'black:','linewidth',1.5)
% plot(C(:,1)*1e-9,C(:,8),'black:','linewidth',1.5)
% plot(C(:,1),C(:,5),'black:','linewidth',1.5)
% plot(C(:,1)*1e-9,C(:,2),'black:','linewidth',1.5)
% plot(C(:,1)*1e-9,C(:,3),'black:','linewidth',1.5)

% % plot(C(:,1),C(:,8),'-.')
% % plot(C(:,1),C(:,9),'-.')


NumTicks = 5;
L = get(gca,'YLim');
% set(gca,'YTick',linspace(0,0.75e-3,NumTicks))
set(gca,'Fontsize', 14)
xlabel(' $\omega$[GeV] ','Interpreter','latex','fontsize',15)
ylabel('$\omega dP/d\omega$','Interpreter','latex','fontsize',15)
lgd = legend('$\uparrow \uparrow \mathbf{e}_L $','$\downarrow \downarrow \mathbf{e}_L $','$\downarrow \uparrow \mathbf{e}_L $','$\uparrow \uparrow \mathbf{e}_R $','$\downarrow \downarrow \mathbf{e}_R $');

set(lgd,'Fontsize',13)

% axis([0 25 1e-7 5e-3])

%box on
%set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[2*21*0.31, 2*12*0.42],'PaperPosition',[0, 0, 2*21*0.31, 2 * 12*0.42],'Position',[0 0 2*21*0.31, 2*12*0.42])

% saveas(f,'fig2.pdf')

% figure(3)




% figure(3)
% plot(A(:,1),A(:,2)+A(:,6))
% hold all
% plot(A(:,1),A(:,3)+A(:,7),'-.')
% plot(A(:,1),A(:,4)+A(:,8))
% plot(A(:,1),A(:,5)+A(:,9),'-.')

% figure(4)
% plot(A(:,1),A(:,2)+A(:,4)+A(:,6)+A(:,8))
% hold all
% plot(A(:,1),A(:,3)+A(:,5)+A(:,7)+A(:,9))

y=A(:,1)/E0;

% figure
% plot(y,pol1./A(:,1))
% hold all
% plot(y,pol2./A(:,1))
% plot(y,(pol1+pol2)./A(:,1),'linewidth',2)
% plot(C(:,1),0.5*(C(:,2)+C(:,3)+C(:,4)+C(:,5)+C(:,6)+C(:,7)+C(:,8)+C(:,9)),'-.','linewidth',2)
% plot(C(:,1),0.5*(C(:,2)+C(:,3)+C(:,4)+C(:,5)),'-.','linewidth',2)
% plot(C(:,1),0.5*(C(:,6)+C(:,7)+C(:,8)+C(:,9)),'-.','linewidth',2)
% legend('Pol1','Pol2','Sum','Volkov sum','Volkov pol 1','Volkov pol 2')



f2=figure
set(f2,'DefaultTextInterpreter','Latex')
set(f2,'DefaultAxesTickLabelInterpreter','Latex')
#set(f2,'DefaultlegendInterpreter','Latex')
% hold all
F0=0;
F2=0;
for n=1:100
  n
    rn=y*(1+eta^2)./((1-y)*n*x);
    cn=1-2*rn;
    sn=2*sqrt(rn.*(1-rn));
    zn=eta/sqrt(1+eta^2)*n*sn;
    f=besselj(n-1,zn).^2+besselj(n+1,zn).^2-2*besselj(n,zn).^2;
    g=4*n^2*besselj(n,zn).^2./zn.^2;
    h=besselj(n-1,zn).^2-besselj(n+1,zn).^2;
    
    F0=F0+((1./(1-y)+1-y).*f-sn.^2./(1+eta^2).*g).*(rn<1);
    F2=F2+(1./(1-y)+1-y).*cn.*h.*(rn<1);
    
   
%     hold all
end
 plot(A(:,1)/E0,-F2./F0);
 hold all
% plot(D(1,:),0.5*(D(2,:)+D(3,:)+D(4,:)+D(5,:)+D(6,:)+D(7,:)+D(8,:)+D(9,:)),'-.')





NumTicks = 5;
L = get(gca,'YLim');
% set(gca,'YTick',linspace(0,0.75e-3,NumTicks))
set(gca,'Fontsize', 14)
xlabel(' $\omega/\varepsilon$ ','Interpreter','latex','fontsize',15)
ylabel('Degree of polarization','Interpreter','latex','fontsize',15)
lgd = legend('Monochromatic','Monochromatic $0.1/\gamma$ collimation','Pulse','Pulse $0.1/\gamma$ collimation');

set(lgd,'Fontsize',13)
rect = [0.56, 0.23, .1, .1];
set(lgd, 'Position', rect)
axis([0 0.8 -1 1])

%box on
%set(f2, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[2*21*0.31, 2*12*0.42],'PaperPosition',[0, 0, 2*21*0.31, 2 * 12*0.42],'Position',[0 0 2*21*0.31, 2*12*0.42])

saveas(f2,'fig3.pdf')

% N=2;
% phi=linspace(-N*pi,N*pi,300);
% eta=1;
% w0=100;
% Estr=m*w0*eta/e;
% field=Estr*sinh(phi)./cosh(phi).^2;
% chi=2*abs(field)*e/m^2*E/m

% ysynch=0;
% for i=1:length(phi)
%     ysynch=ysynch+synchpower(A(:,1)',chi(i),E)*(phi(2)-phi(1))/(2*w0);    
% end
% 
% 
% plot(A(:,1)',ysynch,'g')

% figure(6)
% plot(A(:,1),(pol1-pol2)./(pol1+pol2))
% 
% figure(7)
% plot(B(:,8),B(:,4))