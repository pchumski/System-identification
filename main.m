close all; clear; clc;
 
%% Zauicjowanie zmiennych globalnych
Tp = 1;
 
%% Wczytanie danych pomiarowych
data = readtable('robot_arm.dat');
u = data{1:512, 1};
y = data{1:512, 2};
N = length(u);
t = (0:length(u) - 1) * Tp;
%% Wykresy
figure()
subplot(3,1,1)
plot(t,u)
title('sygnał wejściowy u')
subplot(3,1,2)
plot(t,y)
title('sygnał wyjściowy y')
subplot(3,1,3)
plot(t,u,t,y);
legend('u','y')
%% Metoda korelacyjna
 
M = 126;
Ruu = zeros(M,M);
ryu = zeros(M,1);
tn = 0:1:M-1;
for i=1:M
    for j=1:M
        Ruu(i,j) = Covar([u,u],j-i);
    end
    ryu(i,1) = Covar([y,u],i-1); 
end
gM = inv(Ruu)*ryu;
hM = cumsum(gM * Tp);
figure()
subplot(1,2,1)
plot(tn*Tp,gM);
title('Odpowiedz impulsowa')
subplot(1,2,2)
plot(tn, hM);
title('Odpowiedź skokowa');
 
%% Metoda analizy widmowej
[ruu, tauu] = xcorr(u,'biased');
[ryu, tauuy] = xcorr(y,u,'biased');
Mww = N/10;
Mw = ceil(Mww);
k = 0:1:N-1;
Np = (N/2)-1;
omegak = 2*pi*k/N; % [rad]
omega = omegak/Tp; %[rad/s]
wh = zeros(N,1);
% Okno Hamminga
for i=0:N-1
    if (i < Mw)
        wh(i+1,1) = 0.5*(1 + cos(i*pi/Mw));
    end
end
ru = zeros(N,1);
ry = zeros(N,1);
for i=0:N-1
    ru(i+1,1) = Covar([u,u],i)*wh(i+1);
    ry(i+1,1) = Covar([y,u],i)*wh(i+1);
end
 
Gn = fft(y)./fft(u);
Gw = fft(ryu(N:end))./fft(ruu(N:end));
mod = abs(Gn(1:Np));
arg = angle(Gn(1:Np))*180/pi;
for i = 1:Np
    if arg(i)>2*pi
        arg(i) = arg(i)-360;
    end
end
mod1 = abs(Gw(1:Np));
arg1 = angle(Gw(1:Np))*180/pi;
for i = 1:Np
    if arg1(i)>2*pi
        arg1(i) = arg1(i)-360;
    end
end
Lm = 20*log10(mod);
Lm1 = 20*log10(mod1);
figure()
subplot(1,2,1)
semilogx(omega(1:Np), Lm, 'b')
title('Moduł')
subplot(1,2,2)
semilogx(omega(1:Np),arg,'k')
title('Kat')
%% Podział danych pomiarowych na 2 podzbiory
N_est = fix(N/2);
N_wer = length(u) - N_est;
 
u_est = u(1:N_est);
u_wer = u(N_est+1:end);
 
y_est = y(1:N_est);
y_wer = y(N_est+1:end);
%% Model oscylacyjny
fi = zeros(length(u_est), 4);
for i = 5:length(y_est)
    fi(i,1) = -1*y_est(i-1);
    fi(i,2) = -1*y_est(i-2);
    fi(i,3) = -1*y_est(i-3);
    fi(i,4) = u_est(i-4);
end
p = pinv(fi) * y_est;
ym = zeros(1,length(y_wer));
yp = zeros(1,length(y_wer));
ym(1) = -0.234255; % wartość początkowa symulator
yp(1) = -0.234255; % wartość początkowa predyktor
ym(2) = -0.234255; % wartość początkowa symulator
yp(2) = -0.234255; % wartość początkowa predyktor
ym(3) = -0.234255; % wartość początkowa symulator
yp(3) = -0.234255; % wartość początkowa predyktor
ym(4) = -0.234255; % wartość początkowa symulator
yp(4) = -0.234255; % wartość początkowa predyktor
 
for i = 5:length(y_wer)
    ym(i) = -1*p(1)*ym(i-1) - p(2)*ym(i-2) - p(3)*ym(i-3)  + p(4)*u_wer(i-4);
    yp(i) = -1*p(1)*y_wer(i-1)- p(2)*y_wer(i-2) - p(3)*y_wer(i-3) + p(4)*u_wer(i-4);
end
figure()
hold on
plot(y_wer)
plot(yp,'k')
plot(ym)
title('Model oscylacyjny')
legend('zmierzone','odpowiedz predykowana','odpowiedz symulatora')
hold off
 
%% Identyfikacja - rząd 1
% phiT(n) = [-y(n-1) u(n-1)]
Phi1 = zeros(N_est, 2);
 
for k=2:N_est
    Phi1(k, :) = [-y_est(k-1), u_est(k-1)];
end
 
pLS1 = pinv(Phi1) * y_est;
G1 = tf(pLS1(2), [1, pLS1(1)], Tp);
 
y_hat1 = zeros(N_wer, 1);
 
for k=2:N_wer
    y_hat1(k) = -pLS1(1) * y_wer(k-1) + pLS1(2) * u_wer(k-1);
end
 
%% Identyfikacja - rząd 2
% phiT(n) = [-y(n-1) -y(n-2) u(n-1) u(n-2)]
Phi2 = zeros(N_est, 4);
Phi2(2, :) = [-y_est(1), 0, u_est(1), 0];
 
for k=3:N_est
    Phi2(k, :) = [-y_est(k-1), -y_est(k-2), u_est(k-1), u_est(k-2)];
end
 
pLS2 = pinv(Phi2) * y_est;
G2 = tf([pLS2(3), pLS2(4)], [1, pLS2(1), pLS2(2)], Tp);
 
y_hat2 = zeros(N_wer, 1);
y_hat2(2) = -pLS2(1) * y_wer(1) + pLS2(3) * u_wer(1);
 
for k=3:N_wer
    y_hat2(k) = -pLS2(1) * y_wer(k-1) -pLS2(2) * y_wer(k-2) + pLS2(3) * u_wer(k-1) + pLS2(4) * u_wer(k-2);
end
 
%% Identyfikacja - rząd 3
% phiT(n) = [-y(n-1) -y(n-2) -y(n-3) u(n-1) u(n-2) u(n-3)]
Phi3 = zeros(N_est, 6);
Phi3(2, :) = [-y_est(1), 0, 0, u_est(1), 0, 0];
Phi3(3, :) = [-y_est(2), -y_est(1), 0, u_est(2), u_est(1), 0];
 
for k=4:N_est
    Phi3(k, :) = [-y_est(k-1), -y_est(k-2), -y_est(k-3), u_est(k-1), u_est(k-2), u_est(k-3)];
end
 
pLS3 = pinv(Phi3) * y_est;
G3 = tf(pLS3(4:6)', [1, pLS3(1:3)'], Tp);
 
y_hat3 = zeros(N_wer, 1);
y_hat3(1) = y_wer(1);
y_hat3(2) = -pLS3(1) * y_wer(1) + pLS3(4) * u_wer(1);
y_hat3(3) = -pLS3(1) * y_wer(2) - pLS3(2) * y_wer(1) + pLS3(4) * u_wer(2) + pLS3(5) * u_wer(1);
 
for k=4:N_wer
    y_hat3(k) = -pLS3(1) * y_wer(k-1) - pLS3(2) * y_wer(k-2) - pLS3(3) * y_wer(k-3) + pLS3(4) * u_wer(k-1) + pLS3(5) * u_wer(k-2) + pLS3(6) * u_wer(k-3);
end
 
%% Identyfikacja - rząd 4
% phiT(n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) u(n-1) u(n-2) u(n-3) u(n-4)]
Phi4 = zeros(N_est, 8);
Phi4(2, :) = [-y_est(1), 0, 0, 0, u_est(1), 0, 0, 0];
Phi4(3, :) = [-y_est(2), -y_est(1), 0, 0, u_est(2), u_est(1), 0, 0];
Phi4(4, :) = [-y_est(3), -y_est(2), -y_est(1), 0, u_est(3), u_est(2), u_est(1), 0];
 
for k=5:N_est
    Phi4(k, :) = [-y_est(k-1), -y_est(k-2), -y_est(k-3), -y_est(k-4), u_est(k-1), u_est(k-2), u_est(k-3), u_est(k-4)];
end
 
pLS4 = pinv(Phi4) * y_est;
G4 = tf(pLS4(5:8)', [1, pLS4(1:4)'], Tp);
 
y_hat4 = zeros(N_wer, 1);
y_hat4(1) = y_wer(1);
y_hat4(2) = -pLS4(1) * y_wer(1) + pLS4(5) * u_wer(1);
y_hat4(3) = -pLS4(1) * y_wer(2) - pLS4(2) * y_wer(1) + pLS4(5) * u_wer(2) + pLS4(6) * u_wer(1);
y_hat4(4) = -pLS4(1) * y_wer(3) - pLS4(2) * y_wer(2) - pLS4(3) * y_wer(1) + pLS4(5) * u_wer(3) + pLS4(6) * u_wer(2) + pLS4(7) * u_wer(1);
 
for k=5:N_wer
    y_hat4(k) = -pLS4(1) * y_wer(k-1) - pLS4(2) * y_wer(k-2) - pLS4(3) * y_wer(k-3) - pLS4(4) * y_wer(k-4) + pLS4(5) * u_wer(k-1) + pLS4(6) * u_wer(k-2) + pLS4(7) * u_wer(k-3) + pLS4(8) * u_wer(k-4);
end
 
%% Identyfikacja - rząd 5
% phiT(n) = [-y(n-1) -y(n-2) -y(n-3) -y(n-4) -y(n-5) u(n-1) u(n-2) u(n-3) u(n-4) u(n-5)]
Phi5 = zeros(N_est, 10);
Phi5(2, :) = [-y_est(1), 0, 0, 0, 0, u_est(1), 0, 0, 0, 0];
Phi5(3, :) = [-y_est(2), -y_est(1), 0, 0, 0, u_est(2), u_est(1), 0, 0, 0];
Phi5(4, :) = [-y_est(3), -y_est(2), -y_est(1), 0, 0, u_est(3), u_est(2), u_est(1), 0, 0];
Phi5(5, :) = [-y_est(4), -y_est(3), -y_est(2), -y_est(1), 0, u_est(4), u_est(3), u_est(2), u_est(1), 0];
 
for k=6:N_est
    Phi5(k, :) = [-y_est(k-1), -y_est(k-2), -y_est(k-3), -y_est(k-4), -y_est(k-5), u_est(k-1), u_est(k-2), u_est(k-3), u_est(k-4), u_est(k-5)];
end
 
pLS5 = pinv(Phi5) * y_est;
G5 = tf(pLS5(6:10)', [1, pLS5(1:5)'], Tp);
 
y_hat5 = zeros(N_wer, 1);
y_hat5(1) = y_wer(1);
y_hat5(2) = -pLS5(1) * y_wer(1) + pLS5(6) * u_wer(1);
y_hat5(3) = -pLS5(1) * y_wer(2) - pLS5(2) * y_wer(1) + pLS5(6) * u_wer(2) + pLS5(7) * u_wer(1);
y_hat5(4) = -pLS5(1) * y_wer(3) - pLS5(2) * y_wer(2) - pLS5(3) * y_wer(1) + pLS5(6) * u_wer(3) + pLS5(7) * u_wer(2) + pLS5(8) * u_wer(1);
y_hat5(5) = -pLS5(1) * y_wer(4) - pLS5(2) * y_wer(3) - pLS5(3) * y_wer(2) - pLS5(4) * y_wer(1) + pLS5(6) * u_wer(4) + pLS5(7) * u_wer(3) + pLS5(8) * u_wer(2) + pLS5(9) * u_wer(1);
 
for k=6:N_wer
    y_hat5(k) = -pLS5(1) * y_wer(k-1) - pLS5(2) * y_wer(k-2) - pLS5(3) * y_wer(k-3) - pLS5(4) * y_wer(k-4) - pLS5(5) * y_wer(k-5) + pLS5(6) * u_wer(k-1) + pLS5(7) * u_wer(k-2) + pLS5(8) * u_wer(k-3) + pLS5(9) * u_wer(k-4) + pLS5(10) * u_wer(k-5);
end
 
%% Porównanie odpowiedzi
t_wer = (0:N_wer-1)*Tp;
y1 = lsim(G1, u_wer, t_wer);
y2 = lsim(G2, u_wer, t_wer);
y3 = lsim(G3, u_wer, t_wer);
y4 = lsim(G4, u_wer, t_wer);
y5 = lsim(G5, u_wer, t_wer);
 
% rzad 1
figure()
hold on
plot(y_wer)
plot(y1,'k')
plot(y_hat1)
title('Model ARX 2 parametry LS 1 rzad')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
%rzad 2
figure()
hold on
plot(y_wer)
plot(y2,'k')
plot(y_hat2)
title('Model ARX 4 parametry LS 2 rzad')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
%rzad 3
figure()
hold on
plot(y_wer)
plot(y3,'k')
plot(y_hat3)
title('Model ARX 6 parametry LS 3 rzad')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
%rzad 4
figure()
hold on
plot(y_wer)
plot(y4,'k')
plot(y_hat4)
title('Model ARX 8 parametry LS 4 rzad')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
%rzad5
figure()
hold on
plot(y_wer)
plot(y5,'k')
plot(y_hat5)
title('Model ARX 10 parametry LS 5 rzad')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
 
%% Zmienne instrumentalne
% yest = y(1:256);
% uest = u(1:256);
% ywer = y(257:512);
% uwer = u(257:512);
 
uest = u(1:N_est);
uwer = u(N_est+1:end);
yest = y(1:N_est);
ywer = y(N_est+1:end);
 
x = zeros(1,length(yest));
z = zeros(length(yest),2);
% IV 2 parametry
for i=2:length(yest)
    x(i) = -1*pLS1(1)*x(i-1)+pLS1(2)*uest(i-1);
    z(i,1) = -1*x(i-1);
    z(i,2) = uest(i-1); 
end
piv = (inv(z' * Phi1))*z' * yest;
ymiv = zeros(1,length(ywer));
ypiv = zeros(1,length(ywer));
ymiv(1) = -0.234255;
ypiv(1) = -0.234255;
for i = 2:length(ywer)
    ymiv(i) = -1*piv(1)*ymiv(i-1) + piv(2)*uwer(i-1);
    ypiv(i) = -1*piv(1)*ywer(i-1) + piv(2)*uwer(i-1);
end
figure()
hold on
plot(ywer)
plot(ymiv,'k')
plot(ypiv)
title('IV 2 parametry')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
% IV 4 parametry
x2 = zeros(1,length(yest));
z2 = zeros(length(yest),4);
for i=3:length(yest)
    x2(i) = -1*pLS2(1)*x2(i-1) - pLS2(2)*x2(i-2)+ pLS2(3)*uest(i-1) + pLS2(4)*uest(i-2);
    z2(i,1) = -1*x2(i-1);
    z2(i,2) = -1*x2(i-2);
    z2(i,3) = uest(i-1); 
    z2(i,4) = uest(i-2); 
end
piv2 = (inv(z2' * Phi2))*z2' * yest;
ymiv2 = zeros(1,length(ywer));
ypiv2 = zeros(1,length(ywer));
ymiv2(1) = -0.234255;
ypiv2(1) = -0.234255;
ymiv2(2) = -0.234255;
ypiv2(2) = -0.234255;
for i = 3:length(ywer)
    ymiv2(i) = -1*piv2(1)*ymiv2(i-1) - piv2(2)*ymiv2(i-2) + piv2(3)*uwer(i-1) + piv2(4)*uwer(i-2);
    ypiv2(i) = -1*piv2(1)*ywer(i-1) - piv2(2)*ywer(i-2)  + piv2(3)*uwer(i-1) + piv2(4)*uwer(i-2);
end
figure()
hold on
plot(ywer)
plot(ymiv2,'k')
plot(ypiv2)
title('IV 4 parametry')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
% IV 6 parametry
x3 = zeros(1,length(yest));
z3 = zeros(length(yest),6);
for i=4:length(yest)
    x3(i) = -1*pLS3(1)*x3(i-1) - pLS3(2)*x3(i-2) - pLS3(3)*x3(i-3) + pLS3(4)*uest(i-1) + pLS3(5)*uest(i-2) + pLS3(6)*uest(i-3);
    z3(i,1) = -1*x3(i-1);
    z3(i,2) = -1*x3(i-2);
    z3(i,3) = -1*x3(i-3);
    z3(i,4) = uest(i-1); 
    z3(i,5) = uest(i-2); 
    z3(i,6) = uest(i-3); 
end
piv3 = (inv(z3' * Phi3))*z3' * yest;
ymiv3 = zeros(1,length(ywer));
ypiv3 = zeros(1,length(ywer));
ymiv3(1) = -0.234255;
ypiv3(1) = -0.234255;
ymiv3(2) = -0.234255;
ypiv3(2) = -0.234255;
ymiv3(3) = -0.234255;
ypiv3(3) = -0.234255;
for i = 4:length(ywer)
    ymiv3(i) = -1*piv3(1)*ymiv3(i-1) - piv3(2)*ymiv3(i-2) - piv3(3)*ymiv3(i-3) + piv3(4)*uwer(i-1) + piv3(5)*uwer(i-2) + piv3(6)*uwer(i-3);
    ypiv3(i) = -1*piv3(1)*ywer(i-1) - piv3(2)*ywer(i-2)- piv3(3)*ywer(i-3)  + piv3(4)*uwer(i-1) + piv3(5)*uwer(i-2)+ piv3(6)*uwer(i-3);
end
figure()
hold on
plot(ywer)
plot(ymiv3,'k')
plot(ypiv3)
title('IV 6 parametry')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
 
% IV 8 parametry
x4 = zeros(1,length(yest));
z4 = zeros(length(yest),8);
for i=6:length(yest)
    x4(i) = -1*pLS4(1)*x4(i-1) - pLS4(2)*x4(i-2) - pLS4(3)*x4(i-3) - pLS4(4)*x4(i-4) + pLS4(5)*uest(i-1) + pLS4(6)*uest(i-2) + pLS4(7)*uest(i-3) + pLS4(8)*uest(i-4);
    z4(i,1) = -1*x3(i-1);
    z4(i,2) = -1*x3(i-2);
    z4(i,3) = -1*x3(i-3);
    z4(i,4) = -1*x3(i-4);
    z4(i,5) = uest(i-1); 
    z4(i,6) = uest(i-2); 
    z4(i,7) = uest(i-3); 
    z4(i,8) = uest(i-4); 
end
piv4 = (inv(z4' * Phi4))*z4' * yest;
ymiv4 = zeros(1,length(ywer));
ypiv4 = zeros(1,length(ywer));
ymiv4(1) = -0.234255;
ypiv4(1) = -0.234255;
ymiv4(2) = -0.234255;
ypiv4(2) = -0.234255;
ymiv4(3) = -0.234255;
ypiv4(3) = -0.234255;
ymiv4(4) = -0.234255;
ypiv4(4) = -0.234255;
ymiv4(5) = -0.234255;
ypiv4(5) = -0.234255;
for i = 6:length(ywer)
    ymiv4(i) = -1*piv4(1)*ymiv4(i-1) - piv4(2)*ymiv4(i-2) - piv4(3)*ymiv4(i-3) - piv4(4)*ymiv4(i-4) + piv4(5)*uwer(i-1) + piv4(6)*uwer(i-2) + piv4(7)*uwer(i-3) + piv4(8)*uwer(i-4);
    ypiv4(i) = -1*piv4(1)*ywer(i-1) - piv4(2)*ywer(i-2) - piv4(3)*ywer(i-3) - piv4(4)*ywer(i-4) + piv4(5)*uwer(i-1) + piv4(6)*uwer(i-2) + piv4(7)*uwer(i-3) + piv4(8)*uwer(i-4);
end
figure()
hold on
plot(ywer)
plot(ymiv4,'k')
plot(ypiv4)
title('IV 8 parametry')
legend('zmierzone','odpowiedz symulatora','odpowiedz predyktora')
hold off
%% pomiar bledow
%średnia wartość kwadratu błędu wyjściowego
oscylacyjny = 0; S_LS_2 = 0; S_LS_4 = 0; S_LS_6 = 0; S_LS_8 = 0; S_IV_2 = 0; S_IV_4 = 0; S_IV_6 = 0;
 
for i = 1:length(y_wer)
    oscylacyjny = oscylacyjny + (y_wer(i) - yp(i))^2;
    S_LS_2 = S_LS_2 + (y_wer(i) - y1(i))^2;
    S_LS_4 = S_LS_4 + (y_wer(i) - y2(i))^2;
    S_LS_6 = S_LS_6 + (y_wer(i) - y3(i))^2;
    S_LS_8 = S_LS_8 + (y_wer(i) - y4(i))^2;
    S_IV_2 = S_IV_2 + (y_wer(i) - ymiv(i))^2;
    S_IV_4 = S_IV_4 + (y_wer(i) - ymiv2(i))^2;
    S_IV_6 = S_IV_6 + (y_wer(i) - ymiv3(i))^2;
end
oscylacyjny = oscylacyjny / length(y_wer);
S_LS_2 = S_LS_2 / length(y_wer);
S_LS_4 = S_LS_4 / length(y_wer);
S_LS_6 = S_LS_6 / length(y_wer);
S_LS_8 = S_LS_8 / length(y_wer);
S_IV_2 = S_IV_2 / length(y_wer);
S_IV_4 = S_IV_4 / length(y_wer);
S_IV_6 = S_IV_6 / length(y_wer);
 
%ls i oscylacyjny bledy estymacji
Vp = zeros(1,5);
Vm = zeros(1,5);
 
ep1 = y_wer' - yp;
Vp(1) = (1/length(y_wer)) * ep1*ep1';
em1 = y_wer' - ym;
Vm(1) = (1/length(y_wer)) * em1*em1';
 
% ep = ywer' - yp2;
% Vp(1) = (1/length(ywer)) * ep*ep';
% em = ywer' - ym2;
% Vm(1) = (1/length(ywer)) * em*em';
 
ep2 = y_wer' - y_hat1';
Vp(2) = (1/length(y_wer)) * ep2*ep2';
em2 = y_wer' - y1';
Vm(2) = (1/length(y_wer)) * em2*em2';
 
ep3 = y_wer' - y_hat2';
Vp(3) = (1/length(y_wer)) * ep3*ep3';
em3 = y_wer' - y2';
Vm(3) = (1/length(y_wer)) * em3*em3';
 
ep4 = y_wer' - y_hat3';
Vp(4) = (1/length(y_wer)) * ep4*ep4';
em4 = y_wer' - y3';
Vm(4) = (1/length(y_wer)) * em4*em4';
 
ep5 = y_wer' - y_hat4';
Vp(5) = (1/length(y_wer)) * ep5*ep5';
em5 = y_wer' - y4';
Vm(5) = (1/length(y_wer)) * em5*em5';
 
tm = [0 2 4 6 8];
figure()
plot(tm,log(Vp),tm,log(Vm))
legend('predyktor','symulator');
xlabel('m')
ylabel('ln')
 
%iv bledy estymacji
Vpiv = zeros(1,3);
Vmiv = zeros(1,3);
epiv1 = ywer' - ypiv;
Vpiv(1) = (1/length(ywer)) * epiv1*epiv1';
emiv1 = ywer' - ymiv;
Vmiv(1) = (1/length(ywer)) * emiv1*emiv1';
epiv = ywer' - ypiv2;
Vpiv(2) = (1/length(ywer)) * epiv*epiv';
emiv = ywer' - ymiv2;
Vmiv(2) = (1/length(ywer)) * emiv*emiv';
epiv2 = ywer' - ypiv3;
Vpiv(3) = (1/length(ywer)) * epiv2*epiv2';
emiv2 = ywer' - ymiv3;
Vmiv(3) = (1/length(ywer)) * emiv2*emiv2';
tmiv = [2 4 6];
figure()
plot(tmiv,log(Vpiv),tmiv,log(Vmiv))
legend('predyktor','symulator');
xlabel('m')
ylabel('ln')
 
%% analiza i weryfikacja wybranego modelu
%Macierz kowariancji
eps6 = zeros(1,N_est);
 
for i = 1:N_est
    eps6(i) = y_wer(i) - Phi3(i,:)*pLS3;
end
 
sum_eps_kw6 = 0;
 
for i = 1:N_est
    sum_eps_kw6 = sum_eps_kw6 + eps6(i)^2;
end
 
var_pu6 = 1 / (N_est - 6) * sum_eps_kw6;
 
%regresor stochastyczny
 
sum_rg6 = 0;
 
for i=1:N_est
    sum_rg6 = (Phi3(i,:)'*Phi3(i,:))/(N_est) + sum_rg6;
end
 
Cov_pLS_6_n = var_pu6*sum_rg6^(-1);
 
PU_95_6np1 = [ pLS3(1)-1.96*sqrt(Cov_pLS_6_n(1,1)/N_est) pLS3(1)+1.96*sqrt(Cov_pLS_6_n(1,1)/N_est)];
PU_95_6np2 = [ pLS3(2)-1.96*sqrt(Cov_pLS_6_n(2,2)/N_est) pLS3(2)+1.96*sqrt(Cov_pLS_6_n(2,2)/N_est)];
PU_95_6np3 = [ pLS3(3)-1.96*sqrt(Cov_pLS_6_n(3,3)/N_est) pLS3(3)+1.96*sqrt(Cov_pLS_6_n(3,3)/N_est)];
PU_95_6np4 = [ pLS3(4)-1.96*sqrt(Cov_pLS_6_n(4,4)/N_est) pLS3(4)+1.96*sqrt(Cov_pLS_6_n(4,4)/N_est)];
PU_95_6np5 = [ pLS3(5)-1.96*sqrt(Cov_pLS_6_n(5,5)/N_est) pLS3(5)+1.96*sqrt(Cov_pLS_6_n(5,5)/N_est)];
PU_95_6np6 = [ pLS3(6)-1.96*sqrt(Cov_pLS_6_n(6,6)/N_est) pLS3(6)+1.96*sqrt(Cov_pLS_6_n(6,6)/N_est)];
 
 
%% J
J1 = mean((y_wer - y_hat1) .^ 2);
J2 = mean((y_wer - y_hat2) .^ 2);
J3 = mean((y_wer - y_hat3) .^ 2);
J4 = mean((y_wer - y_hat4) .^ 2);
J5 = mean((y_wer - y_hat5) .^ 2);
 
J1iv = mean((y_wer - ypiv') .^ 2);
J2iv = mean((y_wer - ypiv2') .^ 2);
J3iv = mean((y_wer - ypiv3') .^ 2);
J4iv = mean((y_wer - ypiv4') .^ 2);
 
 
%% Jfit
Jfit1 = (1 - norm(y_wer - y_hat1) / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100;
Jfit2 = (1 - norm(y_wer - y_hat2) / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100;
Jfit3 = (1 - norm(y_wer - y_hat3) / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100;
Jfit4 = (1 - norm(y_wer - y_hat4) / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100;
Jfit5 = (1 - norm(y_wer - y_hat5) / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100;
 
Jfit1iv = (1 - norm(y_wer - ypiv') / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100
Jfit2iv = (1 - norm(y_wer - ypiv2') / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100
Jfit3iv = (1 - norm(y_wer - ypiv3') / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100
Jfit4iv = (1 - norm(y_wer - ypiv4') / norm(y_wer - mean(y_wer)*ones(size(y_wer)))) * 100
%% FPE
V1 = mean((y_est - Phi1 * pLS1) .^ 2);
FPE1 = V1 * (1 + 1 / N_est) / (1 - 1 / N_est);
 
V2 = mean((y_est - Phi2 * pLS2) .^ 2);
FPE2 = V2 * (1 + 2 / N_est) / (1 - 2 / N_est);
 
V3 = mean((y_est - Phi3 * pLS3) .^ 2);
FPE3 = V3 * (1 + 3 / N_est) / (1 - 3 / N_est);
 
V4 = mean((y_est - Phi4 * pLS4) .^ 2);
FPE4 = V4 * (1 + 4 / N_est) / (1 - 4 / N_est);
 
V5 = mean((y_est - Phi5 * pLS5) .^ 2);
FPE5 = V5 * (1 + 5 / N_est) / (1 - 5 / N_est);
 
 
V1iv = mean((y_est - Phi1 * piv) .^ 2);
FPE1iv = V1iv * (1 + 1 / N_est) / (1 - 1 / N_est);
 
V2iv = mean((y_est - Phi2 * piv2) .^ 2);
FPE2iv = V2iv * (1 + 2 / N_est) / (1 - 2 / N_est);
 
V3iv = mean((y_est - Phi3 * piv3) .^ 2);
FPE3iv = V3iv * (1 + 3 / N_est) / (1 - 3 / N_est);
 
V4iv = mean((y_est - Phi4 * piv4) .^ 2);
FPE4iv = V4iv * (1 + 4 / N_est) / (1 - 4 / N_est);
 
%% AIC
 
AIC1 = N_est * log(V1) + 2 * 1;
AIC2 = N_est * log(V2) + 2 * 2;
AIC3 = N_est * log(V3) + 2 * 3;
AIC4 = N_est * log(V4) + 2 * 4;
AIC5 = N_est * log(V5) + 2 * 5;
 
AIC1iv = N_est * log(V1iv) + 2 * 1;
AIC2iv = N_est * log(V2iv) + 2 * 2;
AIC3iv = N_est * log(V3iv) + 2 * 3;
AIC4iv = N_est * log(V4iv) + 2 * 4;
 
%% Współczynnik uwarunkowania macierzy MI
MI1 = 1/N_est * (Phi1' * Phi1);
cond1 = sqrt(max(eig(MI1' * MI1))) / sqrt(mean(eig(MI1' * MI1)));
 
MI2 = 1/N_est * (Phi2' * Phi2);
cond2 = sqrt(max(eig(MI2' * MI2))) / sqrt(mean(eig(MI2' * MI2)));
 
MI3 = 1/N_est * (Phi3' * Phi3);
cond3 = sqrt(max(eig(MI3' * MI3))) / sqrt(mean(eig(MI3' * MI3)));
 
MI4 = 1/N_est * (Phi4' * Phi4);
cond4 = sqrt(max(eig(MI4' * MI4))) / sqrt(mean(eig(MI4' * MI4)));
 
MI5 = 1/N_est * (Phi5' * Phi5);
cond5 = sqrt(max(eig(MI5' * MI5))) / sqrt(mean(eig(MI5' * MI5)));
 
%% Podsumowanie
table([J1; J2; J3; J4; J5],...
      [Jfit1; Jfit2; Jfit3; Jfit4; Jfit5],...
      [cond1; cond2; cond3; cond4; cond5],...
      [FPE1; FPE2; FPE3; FPE4; FPE5],...
      [AIC1; AIC2; AIC3; AIC4; AIC5],...
      'VariableNames', {'J', 'J_FIT', 'cond', 'FPE', 'AIC'},...
      'RowNames', {'LS 2 param.', 'LS 4 param.', 'LS 6 param.', 'LS 8 param.', 'LS 10 param.'})
  
table([J1iv; J2iv; J3iv; J4iv],...
      [Jfit1iv; Jfit2iv; Jfit3iv; Jfit4iv],...
      [cond1; cond2; cond3; cond4],...
      [FPE1iv; FPE2iv; FPE3iv; FPE4iv],...
      [AIC1iv; AIC2iv; AIC3iv; AIC4iv],...
      'VariableNames', {'J', 'J_FIT', 'cond', 'FPE', 'AIC'},...
      'RowNames', {'IV 2 param.', 'IV 4 param.', 'IV 6 param.', 'IV 8 param.'})