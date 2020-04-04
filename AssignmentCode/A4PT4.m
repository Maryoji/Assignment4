%% PART 4

clear all
close all

R1 = 1;
C = 0.25;
R2 = 2;
L = 0.2;
R3 = 10;
a = 100;
R4 = 0.1;
Ro = 1000;

% V = [V1    V2            V3     V4     V5            i1  iL  i3];
G = [-1/R1  1/R1           0       0      0             1   0   0;
    1/R1 (-1/R1)-(1/R2)    0       0      0             0  -1   0;
      0      0            -1/R3    0      0             0   1   0;
      0      0             0   -1/R4    1/R4            0   0   1;
      0      0             0    1/R4  (-1/R4)-(1/Ro)    0   0   0;
      1      0             0       0      0             0   0   0;
      0      1            -1       0      0             0   0   0;
      0      0             0       1      0             0   0   0]


% V =  [V1  V2  V3   V4  V5 i1  iL i3];
Cm =   [-C   C   0   0   0   0   0  0;
         C  -C   0   0   0   0   0  0;
         0   0   0   0   0   0   0  0;
         0   0   0   0   0   0   0  0;
         0   0   0   0   0   0   0  0;
         0   0   0   0   0   0   0  0;
         0   0   0   0   0   0  -L  0;
         0   0   0   0   0   0   0  0]
 
 
 %B = [0 0 0 0 a*(V3/R3) 0 0 0]
 
  I = 0;
  vin = -10:1:10;  
  V3 = zeros(size(vin));
  
for Vin = -10:10 
    I = I + 1;
    F = [0 0 0 0 0 Vin 0 0];
    B = [0 0 0 0 a*(V3(I)/R3) 0 0 0];
    V = ((G+Cm)\F')+B;
    
    
    V3(I) = V(3);
    Vo(I) = V(5);
end

figure(10)
plot(vin,V3)
title('V3')
figure(11)
plot(vin,Vo)
title('Vo')