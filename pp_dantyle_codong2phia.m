clear all;
clc;

% Wmt = 4*9.8;%gia toc phap tuyen cua muc tieu

Tetamt(1) = 90/57.3;% goc huong ban dau cua muc tieu [rad]
dt = 0.001;
i = 2;

% muc tieu co dong 
t0cd = 1.5;   %thoi diem bat dau co dong
too = 15;   %thoi gian co dong
w = 2*pi/3; %tan so co dong cua muc tieu


%tham so muc tieu
Xmt(1) = 7000;
Ymt(1) = 3000;
Vmt = -300;
% Tetamtc = Wmt/Vmt;
Vxmt(1) = Vmt*cos(Tetamt(1));
Vymt(1) = Vmt*sin(Tetamt(1));
% tham so ten lua
Xtl(1) = 0;
Ytl(1) = 0;
Vtl = 600;
Dx(1) = (Xmt(1) - Xtl(1));
Dy(1) = (Ymt(1) - Ytl(1));
D(1) = sqrt(Dx(1)^2 + Dy(1)^2);
xikma(1) = atan(Dy(1)/Dx(1));                       % goc dg ngam
Qtl(1) = asin(Vmt*sin(Tetamt(1) - xikma(1))/Vtl);   % goc hop(dg ngam,Vtl)
Tetatl(1) = atan((Ymt(1)-Ytl(1))/(Xmt(1)-Xtl(1)));%Qtl(1) +xikma(1);
Vxtl(1) = Vtl*cos(Tetatl(1));
Vytl(1) = Vtl*sin(Tetatl(1));
Vx(1) = Vxmt(1) - Vxtl(1);
Vy(1) = Vymt(1) - Vytl(1);

%-------------------------------------------------------------------

N = 3;    % gia tri ty le tao gia phap tuyen ten lua thuc
Vtc(1) =  -(Dx(1)*Vx(1) + Dy(1)*Vy(1))/D(1);
xikmac(1) = (Dx(1)*Vy(1) - Dy(1)*Vx(1))/(D(1)^2);
Wtlv(1) = 0;
Wtl_(1) = 0;
Wtl(1) = 0;
Tetatlc(1) = Wtl(1)/Vtl;
wy(1) = 0; 
Tetamtc(1) = 30/57.3;
g = 9.8;
nuy_tl(1) = Wtl(1)/g;

%khi dong li tuong
% af1 = 5.134e-05 ;
% af2 = 0.0001027 ;
% af3 = 5.134e-05;
% 
% bf1 =  1;
% bf2 = -1.979;
% bf3 =  0.9791;

% %khi dong 0.001
% af1 = 8.12e-05 ;
% af2 = 0.0001624 ;
% af3 = 8.12e-05;
% 
% bf1 =  1;
% bf2 = -1.983;
% bf3 =  0.9831;

% gaz dong 0.001
af1 = 0.0009263;
af2 = 0.001853;
af3 = 0.0009263;

bf1 =  1;
bf2 = -1.873;
bf3 =  0.8774;

while D(i-1) > 30
    if i*dt >= t0cd && i*dt <= (t0cd+too)
        wy(i) = 10*9.8*cos(w*dt*i);
    else wy(i) = 0;
    end
    if i <=2 
        Tetamtc(i) = wy(i)/Vmt;
    Xmt(i) = Xmt(i-1) + Vxmt(i-1)*dt;
    Ymt(i) = Ymt(i-1) + Vymt(i-1)*dt;
    Tetamt(i) = Tetamt(i-1) + Tetamtc(i)*dt;
    Vxmt(i) = Vmt*sin(Tetamt(i));
    Vymt(i) = Vmt*cos(Tetamt(i));
    % tham so ten lua
    Xtl(i) = Xtl(i-1) + Vxtl(i-1)*dt;
    Ytl(i) = Ytl(i-1) + Vytl(i-1)*dt;
    Tetatl(i) = Tetatl(i-1) + Tetatlc(i-1)*dt;
    Vxtl(i) = Vtl*cos(Tetatl(i));
    Vytl(i) = Vtl*sin(Tetatl(i));
    Vx(i) = Vxmt(i) - Vxtl(i);
    Vy(i) = Vymt(i) - Vytl(i);
    
    Dx(i) = (Xmt(i) - Xtl(i));
    Dy(i) = (Ymt(i) - Ytl(i));
    D(i) = sqrt(Dx(i)^2 + Dy(i)^2);
    xikma(i) = atan(Dy(i)/Dx(i));                       % goc dg ngam
    Qtl(i) = asin(Vmt*sin(Tetamt(i) - xikma(i))/Vtl);   % goc hop(dg ngam,Vtl)
    
    %--------------------------------------------------
    
    Vtc(i) =  -(Dx(i)*Vx(i) + Dy(i)*Vy(i))/D(i);
    xikmac(i) = (Dx(i)*Vy(i) - Dy(i)*Vx(i))/(D(i)^2);
    Wtlv(i) = 0;
    Wtl_(i) = 0;
    Wtl(i) = 0;
    Tetatlc(i) = Wtl(i)/Vtl;
     nuy_tl(i) = Wtl(i)/g;
    else
        Tetamtc(i) = wy(i)/Vmt;
    Xmt(i) = Xmt(i-1) + Vxmt(i-1)*dt;
    Ymt(i) = Ymt(i-1) + Vymt(i-1)*dt;
    Tetamt(i) = Tetamt(i-1) + Tetamtc(i)*dt;
    Vxmt(i) = Vmt*sin(Tetamt(i));
    Vymt(i) = Vmt*cos(Tetamt(i));
    % tham so ten lua
    Xtl(i) = Xtl(i-1) + Vxtl(i-1)*dt;
    Ytl(i) = Ytl(i-1) + Vytl(i-1)*dt;
    Tetatl(i) = Tetatl(i-1) + Tetatlc(i-1)*dt;
    Vxtl(i) = Vtl*cos(Tetatl(i));
    Vytl(i) = Vtl*sin(Tetatl(i));
    Vx(i) = Vxmt(i) - Vxtl(i);
    Vy(i) = Vymt(i) - Vytl(i);
    
    Dx(i) = (Xmt(i) - Xtl(i));
    Dy(i) = (Ymt(i) - Ytl(i));
    D(i) = sqrt(Dx(i)^2 + Dy(i)^2);
    xikma(i) = atan(Dy(i)/Dx(i));                       % goc dg ngam
    Qtl(i) = asin(Vmt*sin(Tetamt(i) - xikma(i))/Vtl);   % goc hop(dg ngam,Vtl)
    
    %--------------------------------------------------
    
    Vtc(i) =  -(Dx(i)*Vx(i) + Dy(i)*Vy(i))/D(i);
    xikmac(i) = (Dx(i)*Vy(i) - Dy(i)*Vx(i))/(D(i)^2);
    Wtlv(i) = N*Vtc(i)*xikmac(i);
%     Wtl(i) = Wtlv(i)/cos(Qtl(i));
    Wtl_(i) = (Wtlv(i)*af1 + Wtlv(i-1)*af2 + Wtlv(i-2)*af3 - Wtl_(i-1)*bf2 - Wtl_(i-2)*bf3)/(bf1);
    Wtl(i) = Wtl_(i)/cos(Qtl(i));
%     if Wtl(i) > 10*9.8
%         Wtl(i) = 10*9.8;
%     end
%     if Wtl(i) < -10*9.8
%         Wtl(i) = -10*9.8;
%     end
    Tetatlc(i) = Wtl(i)/Vtl;
     nuy_tl(i) = Wtl(i)/g;
    end
    
    h(i) = D(i)*D(i)*xikmac(i)/abs(Vtc(i));
    i = i + 1;
   %  nuy_tl(i) = Wtl(i)/g;
end
figure(1)%,'quy dao ten lua muc tieu');
plot(Xmt,Ymt,Xtl,Ytl);
% ylim([0 4000]);
legend('quy dao muc tieu', 'quy dao ten lua')
grid on
hold on

figure(2)
plot(h)
title('Do truot')
grid on
hold on

figure(3)
plot(nuy_tl);
title('Qua tai TL');
legend('QTTL');
grid on;
hold on;