% BKK95_IR.m
% plots impulse responses
% 
% Katrin Rabitsch, March 2008

global T shock;

% x = [B K1 K2 Z1 Z2];
% y = [C1 C2 N1 N2 ARM1 ARM2 X1 X2 P RX A1 B1 A2 B2 W1 W2 R1 R2 LAM1 LAM2 QA1 QB1 QA2 QB2 QQ IR1 GDP1 NX1];  

iC1=1;      iC2=2;
iN1=3;      iN2=4;
iARM1=5;    iARM2=6;
iX1=7;      iX2=8;
iP=9;       iRX=10;
iA1=11;     iB1=12;
iA2=13;     iB2=14;
iW1=15;     iW2=16;
iR1=17;     iR2=18;

iNX1=size(gx,1);     iGDP=size(gx,1)-1;
iIR1=size(gx,1)-2;

% lag var
iB=size(gx,1)+size(hx,1)-6;
iK1=size(gx,1)+size(hx,1)-5;
iK2=size(gx,1)+size(hx,1)-4;
iG1=size(gx,1)+size(hx,1)-3;
iG2=size(gx,1)+size(hx,1)-2;
iZ1=size(gx,1)+size(hx,1)-1;
iZ2=size(gx,1)+size(hx,1);


if shock==1;     FigNameStrg='Impulse response to 1% government shock in country 1'; 
elseif shock==2; FigNameStrg='Impulse response to 1% government shock in country 2'; 
elseif shock==3; FigNameStrg='Impulse response to 1% productivity shock in country 1'; 
elseif shock==4; FigNameStrg='Impulse response to 1% productivity shock in country 2'; 
end



kk=[1:T];
L = 4;   % number of lines of subplots
Col = 2; % number of columns of subplots

figure('Name',FigNameStrg)
subplot(L,Col,1)
plot(kk,IR(1:T,iC1),'b -') 
title('consumption in country 1, C1'),

subplot(L,Col,2)
plot(kk,IR(1:T,iC2),'b -') 
title('consumption in country 2, C2'),

subplot(L,Col,3)
plot(kk,IR(1:T,iX1),'b -') 
title('investment in country 1, X1'),

subplot(L,Col,4)
plot(kk,IR(1:T,iX2),'b -') 
title('investment in country 2, X2'),

subplot(L,Col,5)
plot(kk,IR(1:T,iN1),'b -') 
title('hours worked in country 1, N1'),

subplot(L,Col,6)
plot(kk,IR(1:T,iN2),'b -') 
title('hours worked in country 3, N2'),

subplot(L,Col,7)
plot(kk,IR(1:T,iARM1),'b -') 
title('final output in country 1, ARM1'),

subplot(L,Col,8)
plot(kk,IR(1:T,iARM2),'b -') 
title('final output in country 3, ARM2'),



L = 3;   % number of lines of subplots
Col = 2; % number of columns of subplots

figure('Name',FigNameStrg)
subplot(L,Col,1)
plot(kk,IR(1:T,iK1),'b -') 
title('capital stock in country 1, K1'),

subplot(L,Col,2)
plot(kk,IR(1:T,iK2),'b -') 
title('capital stock in country 2, K2'),

subplot(L,Col,3)
plot(kk,IR(1:T,iW1),'b -') 
title('real wage in country 1, W1'),

subplot(L,Col,4)
plot(kk,IR(1:T,iW2),'b -') 
title('real wage in country 2, W2'),

subplot(L,Col,5)
plot(kk,IR(1:T,iR1),'b -') 
title('real rental rate in country 1, R1'),

subplot(L,Col,6)
plot(kk,IR(1:T,iR2),'b -') 
title('real rental rate in country 2, R2'),




L = 2;   % number of lines of subplots
Col = 2; % number of columns of subplots

figure('Name',FigNameStrg)
subplot(L,Col,1)
plot(kk,IR(1:T,iP),'b -') 
title('Terms of trade, P'),

subplot(L,Col,2)
plot(kk,IR(1:T,iRX),'b -') 
title('Real Exchange Rate, RX'),

subplot(L,Col,3)
plot(kk,IR(1:T,iB),'b -') 
title('Bond holding in country 1, B'),

subplot(L,Col,4)
plot(kk,IR(1:T,iNX1),'b -') 
title('Net Exports in country 1, NX1'),

% subplot(L,Col,5)
% plot(kk,IR(1:T,iA2),'b -') 
% title('exports, A2'),
% 
% subplot(L,Col,6)
% plot(kk,IR(1:T,iB1),'b -') 
% title('imports, B1'),



