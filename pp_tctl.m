clear all;
clc;
%%--------------------------------

g = 9.8;
Wmt = 3*g; % gia toc muc tieu
Tetamt(1) = 180/57.3;
dt  = 0.1;
i = 2;
%% tham so muc tieu
Xmt(1) = 15000;
Ymt(1) = 10000;
Vmt = 600;  %[m/s]
Tetamtc = Wmt/Vmt;
Vxmt(1) = Vmt*cos(Tetamt(1));
Vymt(1) = Vmt*sin(Tetamt(1));
%% tham so ten lua
Xtl(1) = 0;
Ytl(1) = 8000;
Vtl = 1200;
%% tinh toan
Dx(1) = (Xmt(1) - Xtl(1));
Dy(1) = (Ymt(1) - Ytl(1));
D(1) = sqrt(Dx(1)^2 + Dy(1)^2);
zima(1) = atan(Dy(1)/Dx(1)); % goc duong ngam [rad]
Qtl(1) = asin(Vmt*sin(Tetamt(1) - zima(1))/Vtl);    %goc cua (duong ngam, Vtl)
Tetatl(1) = Qtl(1) + zima(1);
Vxtl(1) = Vtl*cos(Tetatl(1));
Vytl(1) = Vtl*sin(Tetatl(1));
Vx(1) = Vxmt(1) - Vxtl(1);
Vy(1) = Vymt(1) - Vytl(1);
%----------------------------------
N = 4; % gia tri ty le tao gia toc phap tuyen ten lua thuc
Vtc(1) = -(Dx(1)*Vx(1) + Dy(1)*Vy(1))/abs(D(1));
zimac(1) = (Dx(1)*Vy(1) - Dy(1)*Vx(1))/(D(1)^2);
Wtlv(1) = N*Vtc(1)*zimac(1);
Wtl(1) = Wtlv(1)/cos(Qtl(1));
Tetatlc(1) = Wtl(1)/Vtl;
while Xmt(i-1) > Xtl(i-1)
    Xmt(i) = Xmt(i-1) + Vxmt(i-1)*dt;
    Ymt(i) = Ymt(i-1) + Vymt(i-1)*dt;
    Tetamt(i) = Tetamt(i-1) + Tetamtc*dt;
    Vxmt(i) = Vmt*cos(Tetamt(i));
    Vymt(i) = -Vmt*sin(Tetamt(i));
    %% tham so ten lua
    Xtl(i) = Xtl(i-1) + Vxtl(i-1)*dt;
    Ytl(i) = Ytl(i-1) + Vytl(i-1)*dt;
    Tetatl(i) = Tetatl(i-1) + Tetatlc(i-1)*dt;
    Vxtl(i) = Vtl*cos(Tetatl(i));
    Vytl(i) = Vtl*sin(Tetatl(i));
    Vx(i) = Vxmt(i) - Vxtl(i);
    Vy(i) = Vymt(i) - Vytl(i);
    
    Dx(i) = Xmt(i) - Xtl(i);
    Dy(i) = Ymt(i) - Ytl(i);
    D(i) = sqrt(Dx(i)^2 + Dy(i)^2);
    zima(i) = atan(Dy(i)/Dx(i));
    Qtl(i) = asin(Vmt*sin(Tetamt(i) - zima(i))/Vtl); 
    %----------------------------------
   
    Vtc(i) = -(Dx(i)*Vx(i) + Dy(i)*Vy(i))/abs(D(i));
    zimac(i) = (Dx(i)*Vy(i) - Dy(i)*Vx(i))/(D(i)^2);
    Wtlv(i) = N*Vtc(i)*zimac(i);
    Wtl(i) = Wtlv(i);%/cos(Qtl(i));
    Tetatlc(i) = Wtl(i)/Vtl;
    i = i + 1;

end
% plot(Tetamt);
% legend('Tetamt');
% hold on
figure();
plot(Xmt,Ymt, Xtl, Ytl);
hold on
grid on
