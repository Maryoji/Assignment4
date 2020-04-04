%%
% *MARY OJI 101036761*
%%
% *CIRCUIT MODELING*
%% PART 1

clear all
close all
%Resistors
R1 = 1;
R2 = 2;
R4 = 0.1;
R3 = 10;
R0 = 1000;

%Capacitor, Inductor and others
C = 0.25;
L = 0.2;
a = 100;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G0 = 1/R0;

G = [ -G1   G1        0     0       0        1   0   0; 
      G1 (-G1)-(G2)   0     0       0        0  -1   0;
      0      0        G3    0       0        0   1   0;
      0      0        0     -G4     G4       0   0   1;  
      0      0        0     G4   (-G4)-(G0)  0   0   0; 
      1      0        0     0       0        0   0   0;
      0      1        -1    0       0        0   0   0; 
      0      0        a/R3  1       0        0   0   0]
  
  
CM =   [-C C 0  0  0  0  0  0;
        C -C 0  0  0  0  0  0;
        0  0 0  0  0  0  0  0; 
        0  0 0  0  0  0  0  0;
        0  0 0  0  0  0  0  0; 
        0  0 0  0  0  0  0  0; 
        0  0 0  0  0  0  -L 0; 
        0  0 0  0  0  0  0  0]
 
i = 0;
for Vin = -10:10 
    i = i + 1;
    F = [0 0 0 0 0 Vin 0 0];
    V = G\F';
    V3(i) = V(3);
    V0(i) = V(5);
end

vin = -10:1:10;
figure(1)
hold on;
title('DC Sweep');
xlabel('Vin sweep from -10 to 10 (V)');
ylabel('Voltage (V)');
plot(vin,V3)
hold on;
plot(vin,V0)
legend('node V3','node V0');

Voi = zeros(1, 1000);
Adb = zeros(1, 1000);
w = 1:1:1000;
for i = 1:1000
    F = [0 0 0 0 0 1 0 0];
    V = (G+(1i*w(i)*CM))\F';
    Voi(i) = V(5);
    Adb(i) = 20*log10(Voi(i)/10);
end

v1 = 1;
figure(2)
semilogx(w,Voi)
title('Vo vs w')
grid on
figure(3)
semilogx(w,Adb)
title('Gain in dB')
grid on

Cnd = 0.25 + 0.05*randn(1,100);
for i = 1:100
    F = [0 0 0 0 0 1 0 0];
    CMND =   [ -Cnd(i)  Cnd(i) 0  0  0  0  0  0;
                Cnd(i) -Cnd(i) 0  0  0  0  0  0;
                0  0 0  0  0  0  0  0; 
                0  0 0  0  0  0  0  0;
                0  0 0  0  0  0  0  0; 
                0  0 0  0  0  0  0  0; 
                0  0 0  0  0  0  L  0; 
                0  0 0  0  0  0  0  0];
    V = (G+(pi*CMND))\F';
    Vh(i) = V(5);
end
figure(4)
hist(Vh)
title('Histogram of Gain')

%% PART 2

VPulse = zeros(8,1);
VSin = zeros(8,1);
VGauss = zeros(8,1);
tstep = 0.001;
time =0;
VinPulse = 0;
deltat = 0.001;
w2 = 2*pi*(1/0.03);
figure(5);
clf;
figure(6);
clf;
J = 0;
for J = 1:1000 %each step represents a milisecond

    if time == 0.03
        VinPulse = 1;
    end
    
   VinSin = sin(w2*time);

   VinGauss = exp(-(time-0.06).^2/(2*(0.03)^2));
    
    FPulse = [0 0 0 0 0 VinPulse 0 0];
    FSin = [0 0 0 0 0 VinSin 0 0];
    FGauss = [0 0 0 0 0 VinGauss 0 0];
    A = G+(CM/deltat);
    
    VPulse = A\(CM*(VPulse/deltat)+FPulse.');
    VSin = A\(CM*(VSin/deltat)+FSin.');
    VGauss = A\(CM*(VGauss/deltat)+FGauss.');
       
    time = tstep*J;
    
    VinPulseIn(J,1) = VinPulse;
    VinSinIn(J,1) = VinSin;
    VinGaussIn(J,1) = VinGauss;
    
    VPulseO(J,1) = VPulse(5);
    VSinO(J,1) = VSin(5);
    VGaussO(J,1) = VGauss(5);
    
end

figure(5)
subplot(3,1,1)
title('Blue-input Red-output')
plot(deltat:deltat:time,VinPulseIn,deltat:deltat:time,VPulseO)
subplot(3,1,2)
plot(deltat:deltat:time,VinSinIn,deltat:deltat:time,VSinO)
subplot(3,1,3)
plot(deltat:deltat:time,VinGaussIn,deltat:deltat:time,VGaussO)


figure(6)
XPulse = abs(fft(VinPulseIn));%fft(Vin11,length(Vin11));
XSin = abs(fft(VinSinIn));
XGauss = abs(fft(VinGaussIn));
XPulseOut = abs(fft(VPulseO));
XSinOut = abs(fft(VSinO));
XGaussOut = abs(fft(VGaussO));
subplot(3,1,1)
plot(-(time/2-deltat):deltat:time/2,XPulse,-(time/2-deltat):deltat:time/2,XPulseOut)%plot(1./(deltat:deltat:time),X)
title('fft Blue-input Red-output')
grid on
subplot(3,1,2)
plot(-(time/2-deltat):deltat:time/2,XSin,-(time/2-deltat):deltat:time/2,XSinOut)
grid on
subplot(3,1,3)
plot(-(time/2-deltat):deltat:time/2,XGauss,-(time/2-deltat):deltat:time/2,XGaussOut)
grid on


figure(7)
XshiftPulse = fftshift(XPulse);
XshiftSin = fftshift(XSin);
XshiftGauss = fftshift(XGauss);
XshiftPulseOut = fftshift(XPulseOut);
XshiftSinOut = fftshift(XSinOut);
XshiftGaussOut = fftshift(XGaussOut);
subplot(3,1,1)
grid on
title('fftshift Blue-input Red-output')
plot(-(time/2-deltat):deltat:time/2,XshiftPulse,-(time/2-deltat):deltat:time/2,XshiftPulseOut)
subplot(3,1,2)
grid on
plot(-(time/2-deltat):deltat:time/2,XshiftSin,-(time/2-deltat):deltat:time/2,XshiftSinOut)
subplot(3,1,3)
plot(-(time/2-deltat):deltat:time/2,XshiftGauss,-(time/2-deltat):deltat:time/2,XshiftGaussOut)
grid on

%%%%%%% As the time step increases the smoother the
%%%%%%% fourier transform plot becomes. Which would mean it is more 
%%%%%%% precises.