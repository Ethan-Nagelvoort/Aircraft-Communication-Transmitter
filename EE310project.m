% Filename:     	EE310project.m
% Author:       	Karl Horcasitas
%               	Ethan Nagelvoort
%               	Giemer Lozares
% Last Modified:	4/18/2020
 
% This script computes the frequency response for a s-domain transfer function.
 
clearvars
close all;
w = (130*10^6*2*pi); %our frequency is 130Mhz * 2pi
 
z1 = 15 + (j*w*24*10^(-9));  %initilize variables to hold our impedance values
z2 = 15 + (j*w*24*10^(-9));
z3 = 1/(j*w*31.5*10^(-12));
z4 = j*w*(24.2*10^(-9));
z5 = j*w*(24.2*10^(-9));
z6 = 1/(j*w*1*10^(-15));
z7 = 1/(j*w*22*10^(-12));
z8 = j*w*(56*10^(-9));
z9 = 1/(j*w*44*10^(-12));
 
%z10:
initL10 = j*w*(56*10^(-9));
initC10 = 1/(j*w*4.7*10^(-12));
z10 = (initL10*initC10)/(initL10+initC10);
 
z11 = 1/(j*w*33*10^(-12));
 
%z12
initL12 = j*w*(39*10^(-9));
initC12 = 1/(j*w*12*10^(-12));
z12 = (initL12*initC12)/(initL12+initC12);
 
z13 = 1/(j*w*16.7*10^(-12));
z14 = 50;
 
v1 = 1; v2 = v1; vR = 0; %using v1=v2=1cos(2pi*130mHz*t), vR = 0
 
A = [(z1+z2+z3) (-z3) 0 0 0 0 0; %Setting up matrix A using our mesh equation
	(-z3) (z3+z4+z5+z6) (-z6) 0 0 0 0;
	0 (-z6) (z6+z7) (-z7) 0 0 0;
	0 0 (-z7) (z7+z8+z9) (-z9) 0 0;
	0 0 0 (-z9) (z9+z10+z11) (-z11) 0;
	0 0 0 0 (-z11) (z11+z12+z13) (-z13);
	0 0 0 0 0 (-z13) (z13+z14)];
B = [(v1+v2);0;0;0;0;0;(-vR)]; %Setting up matrix B (note that vT2-VT1=0 in the third row)
X = linsolve(A,B); %linearly solve the matrices and store in X (X=A^-1*B)
 
%since v2=v1,
%now we raise v1 and v2 until power in z14 is 80 watts.
v1 = 71.180; v2=v1; %using these values of v1,v2, we get a power of 79.9942W
B = [(v1+v2);0;0;0;0;0;(-vR)];
X3 = linsolve(A,B); %using these new currents, we use eq 1/2 * Z14 * abs(i7)^2 to get power
 
ourI5 = X3(5,1); %convert i5 to polar to check answers
rhoI5=abs(ourI5);
thetaI5=angle(ourI5);
%-------------------
ourI7 = X3(7,1); %convert i7 to polar to check answers
rhoI7=abs(ourI7);
thetaI7=angle(ourI7);
 
%----------------------------------------Irms & Vrms Values
Power = 0.5 * (rhoI7^2) * (z14);
 
%RMS VALUES:
%R1:
R1Irms = (abs(X3(1,1)))/(sqrt(2)); %get the real part of i1 and divide by square root 2
R1V = X3(1,1) * 15; %R1 voltage drop = 15 ohm * i1
R1Vrms = abs(R1V)/(sqrt(2)); %get the real part of the R1 voltage drop and divide by square root 2
%R2:
R2Irms = R1Irms; %since R2 = R1, it will have the same Irms
R2Vrms = R1Vrms; %also the same Vrms
%L1:
L1Irms = R1Irms; %since L1 is in series with R1, they will have the same Irms
L1V = X3(1,1) * (j * w * 24*10^-9); %L1 voltage drop is i1 * (jw * 24nano Hz)
L1Vrms = abs(L1V)/(sqrt(2)); %get the real part of VL1 and divide by square root 2
%L2:
L2Irms = L1Irms; %since L2 is equal to L1 and in the same mesh, they have the same Irms
L2Vrms = L1Vrms; %also the same Vrms
%C3:
C3I = (X3(1,1)-X3(2,1)); %capacitor 3 current is (i1-i2)
C3Irms = abs(C3I)/sqrt(2); %get real part of C3I and divide by square root 2
C3V = C3I * (1/(j * w *31.5*10^-12)); %C3 voltage drop is C3I * ZC3
C3Vrms = abs(C3V)/sqrt(2); %get real part of C3 voltage drop and divide by square root 2
%L4:
L4Irms = abs(X3(2,1))/sqrt(2); %get real part of i2 and divide by square root 2
L4Vrms = abs(L4Irms * (j * w * 24.2*10^-9)); %get real part of L4 voltage drop (I * ZL4) and divide by square root 2
%L5:
L5Irms = L4Irms; %since L5 is equal to L4 and in the same mesh, they have the same Irms
L5Vrms = L4Vrms; %also the same Vrms
%C6:
C6Irms = abs(X3(2,1) - X3(3,1))/sqrt(2); %get real part of C6I (i2-i3) and divide by square root 2
C6Vrms = abs(C6Irms*(1/(j * w * 1*10^-15))); %get real part of C6 voltage drop (C6I * ZC6) using Irms
%C7:
C7Irms = abs(X3(3,1) - X3(4,1))/sqrt(2); %get real part of C7I (i3-i4) and divide by square root 2
C7Vrms = abs(C7Irms*(1/(j * w * 22*10^-12))); %get real part of C7 voltage drop (C7I * ZC7) using Irms
%L8:
L8Irms = abs(X3(4,1))/sqrt(2); %get real part of i4 and divide by square root 2
L8Vrms = abs(L8Irms*(j*w*56*10^-9)); %get real part of L8 voltage drop (L8I * ZL8) using Irms
%C9:
C9Irms = abs(X3(4,1) - X3(5,1))/sqrt(2); %get real part of C9I (i4-i5) and divide by square root 2
C9Vrms = abs(C9Irms*(1/(j * w * 44*10^-12))); %get real part of C9 voltage drop (C9I * ZC9) using Irms
%C10:
VZ10 = X3(5,1) * z10; %since we calculate impedance of L10 and C10 in parallel, we must split them up
VZ10rms = abs(VZ10)/sqrt(2); %get real part of that Z10 voltage drop and divide by square root 2
ZC10 = 1/(j*w*4.7*10^-12); %the impedance of C10
C10I = VZ10/ZC10; %since we have VZ10, divide by ZC10 (ohms law) to get current
C10Irms = abs(C10I)/sqrt(2); %get real part of C10I and divide by square root 2
C10Vrms = VZ10rms; %C10 parallel to L10, so they will have the same voltage drop and thus Vrms
%L10:
ZL10 = j * w * 56*10^-9; %impedance of L10
L10I = VZ10/ZL10; %since we have VZ10, divide by ZL10 to get current
L10Irms = abs(L10I)/sqrt(2); %get real part of L10I and divide by square root 2
L10Vrms = VZ10rms; %C10 parallel to L10, so they will have the same voltage drop and thus Vrms
%C11:
C11Irms = abs(X3(5,1) - X3(6,1))/sqrt(2); %get real part of C11 current (i5-i6) and divide by square root 2
C11Vrms = abs(C11Irms*(1/(j * w * 33*10^-12))); %get real part of C11 voltage drop (C11I * ZC11) using Irms
%C12:
VZ12 = X3(6,1) * z12; %like for VZ10, calculate the voltage drop of the parallel components
VZ12rms = abs(VZ12)/sqrt(2); %get real part of VZ12 and divide by square root 2
ZC12 = 1/(j*w*12*10^-12); %the impedance of C12
C12I = VZ12/ZC12; %since we have VZ12, divide by ZC12 (ohms law) to get current
C12Irms = abs(C12I)/sqrt(2); %get real part of C12I and divide by square root 2
C12Vrms = VZ12rms; %C12 parallel to L12, so they will have the same voltage drop and thus Vrms
%L12:
ZL12 = j * w * 39*10^-9; %impedance of L12
L12I = VZ12/ZL12; %since we have VZ12, divide by ZL12 to get current
L12Irms = abs(L12I)/sqrt(2); %get real part of L12I and divide by square root 2
L12Vrms = VZ12rms; %C12 parallel to L12, so they will have the same voltage drop and thus Vrms
%C13:
C13Irms = abs(X3(6,1) - X3(7,1))/sqrt(2); %get real value of C13I (i6-i7) and divide by square root 2
C13Vrms = abs(C13Irms*(1/(j * w * 16.7*10^-12))); %get real part of C13 voltage drop (C13I * ZC13) using Irms
%Z14:
Z14Irms = abs(X3(7,1))/sqrt(2); %get real value of Z14I and divide by square root 2
Z14Vrms = abs(Z14Irms*50); %using ohms law and Irms get Vrms
 
RMSvaluesUnrounded = [R1Vrms R1Irms; %store RMS values in a matrix
	L1Vrms L1Irms;
	R2Vrms R2Irms;
	L2Vrms L2Irms;
	C3Vrms C3Irms;
	L4Vrms L4Irms;
	L5Vrms L5Irms;
	C6Vrms C6Irms;
	C7Vrms C7Irms;
	L8Vrms L8Irms;
	C9Vrms C9Irms;
	L10Vrms L10Irms;
	C10Vrms C10Irms;
	C11Vrms C11Irms;
	L12Vrms L12Irms;
	C12Vrms C12Irms;
	C13Vrms C13Irms;
	Z14Vrms Z14Irms];
 
RMSvalues = round(RMSvaluesUnrounded*100)/100; %round all RMS values to 2 decimal places
 
%------------------------Thevinin Model
v1b=0; v2b=v1b;
%with vR=1(cos(2pi * 130MHz * t), vR=1
vRb=1;
 
Ath = [(z1+z2+z3) (-z3) 0 0 0 0 0; %Matrix for Thevin Voltage
	(-z3) (z3+z4+z5+z6) (-z6) 0 0 0 0;
	0 (-z6) (z6+z7) (-z7) 0 0 0;
	0 0 (-z7) (z7+z8+z9) (-z9) 0 0;
	0 0 0 (-z9) (z9+z10+z11) (-z11) 0;
	0 0 0 0 (-z11) (z11+z12+z13) (-z13);
	0 0 0 0 0 (-z13) (z13+z14)];
Bth = [(v1b+v2b);0;0;0;0;0;(-vRb)];
Xth = linsolve(Ath,Bth);
 
%Vth = i3-i4 * Z7
Vth = (Xth(3,1)-Xth(4,1))*z7;
 
%make C7 so huge so that it behaves like a short
z7th = 1/(j * w * 0.01*10^-6);
Ath2 = [(z1+z2+z3) (-z3) 0 0 0 0 0; %Matrix for Norton Current
	(-z3) (z3+z4+z5+z6) (-z6) 0 0 0 0;
	0 (-z6) (z6+z7th) (-z7th) 0 0 0;
	0 0 (-z7) (z7th+z8+z9) (-z9) 0 0;
	0 0 0 (-z9) (z9+z10+z11) (-z11) 0;
	0 0 0 0 (-z11) (z11+z12+z13) (-z13);
	0 0 0 0 0 (-z13) (z13+z14)];
 
Xth2 = linsolve(Ath2,Bth);
 
Inorton = (Xth2(3,1) - Xth2(4,1)); %Inorton = i3-i4 (the current of C7)
 
Zth = Vth/Inorton; %our Thevinin impedance Vth/Inorton
 
%the conjugate of our thevenin impedance
ZthConj = 19.561400245825563 - 1.015825791894316i;
Rth = 19.561400245825563; %the resistive part of our thevivin impedance
%using eq. 11.20 in textbook
%Pmax = abs(Vth^2)/8Rth
PmaxTh = ((abs(Vth))^2)/(8*Rth);
 
%Complex part of Zth
LforZth = 1.015825791894316i/(i*w); %our thevinin impedance as an inductor
CforZth = 1/(-1.015825791894316i * i * w); %as a capacitor (for our conjugate model)
 
HalfOfVth = Vth/2; %calculate half of Vth to check our work
 
%Vth as sinusoid
VthToSin = -0.311889913008039 + 0.375749839024849i;
rhoVthToSin=abs(VthToSin); %real part of the sinusoid
thetaVthToSin=angle(VthToSin);
VthSin=rhoVthToSin*exp(1i*thetaVthToSin);

