% BKK95_model.m
% 
% March 2008, Katrin Rabitsch

function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f,eta] =  BKK95_model(approx);

global dummy_sig

% approx = 1; %Order of approximation desired

%Define parameters
syms omeg sig betta quc rho_Z1 rho_Z2 delta gam mue1 mue2 alfa phiK rho_Z2Z1 rho_Z1Z2 rho_G1 rho_G2 rho_G1G2 rho_G2G1 G1bar G2bar eta sigmaZ sigmaG

%Define variables 
syms B  K1  K2  Z1  Z2  C1  C2  N1  N2  ARM1  ARM2  X1  X2  P  RX  A1  B1  A2  B2  W1  W2  R1  R2  LAM1  LAM2  QA1  QB1  QA2  QB2  QQ  IR1  GDP1  NX1  G1  G2  K1lead  K2lead
syms Bp K1p K2p Z1p Z2p C1p C2p N1p N2p ARM1p ARM2p X1p X2p Pp RXp A1p B1p A2p B2p W1p W2p R1p R2p LAM1p LAM2p QA1p QB1p QA2p QB2p QQp IR1p GDP1p NX1p G1p G2p K1leadp K2leadp
syms U1 U2 U1p U2p F1 F2 F1p F2p




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% utility and marginal utility
U1  = (1/(1-gam))*(C1 ^mue1*(1-N1 )^(1-mue1))^(1-gam);
U1p = (1/(1-gam))*(C1p^mue1*(1-N1p)^(1-mue1))^(1-gam);
U2  = (1/(1-gam))*(C2 ^mue2*(1-N2 )^(1-mue2))^(1-gam);
U2p = (1/(1-gam))*(C2p^mue2*(1-N2p)^(1-mue2))^(1-gam);

% production function in tradable sector and marginal products
F1  = Z1 *K1 ^alfa*N1 ^(1-alfa);
F1p = Z1p*K1p^alfa*N1p^(1-alfa);
F2  = Z2 *K2 ^alfa*N2 ^(1-alfa);
F2p = Z2p*K2p^alfa*N2p^(1-alfa);

% diff(U1,'C1')
% diff(U2,'C2')
% diff(U1,'N1')
% diff(U2,'N2')
% 
% diff(F1,'N1')
% diff(F1,'K1')
% diff(F2,'N2')
% diff(F2,'K2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write equations, ei=1:

% FINAL GOOD PRODUCERS:

e1 = A1 - (   omeg *(QA1)^(-sig)*ARM1);  % country 1's demand for a good
e2 = B2 - (   omeg *(QB2)^(-sig)*ARM2);  % country 2's demand for b good
e3 = B1 - ((1-omeg)*(QB1)^(-sig)*ARM1);  % country 1's demand for b good
e4 = A2 - ((1-omeg)*(QA2)^(-sig)*ARM2);  % country 2's demand for a good

if dummy_sig==1;
    e5 = ARM1 - (((A1^   omeg) *(B1^(1-omeg)))/((omeg^omeg)*((1-omeg)^(1-omeg)))); % Armington aggregator, final output in country 1
    e6 = ARM2 - (((A2^(1-omeg))*(B2^   omeg) )/((omeg^omeg)*((1-omeg)^(1-omeg)))); % Armington aggregator, final output in country 2  
else;
    e5 = ARM1 - (((   omeg ^(1/sig))*(A1^((sig-1)/sig))+((1-omeg)^(1/sig))*(B1^((sig-1)/sig)))^(sig/(sig-1)));  % Armington aggregator, final output in country 1
    e6 = ARM2 - ((((1-omeg)^(1/sig))*(A2^((sig-1)/sig))+(   omeg ^(1/sig))*(B2^((sig-1)/sig)))^(sig/(sig-1)));  % Armington aggregator, final output in country 2               
end;

% INTERMEDIATE GOOD PRODUCERS:

e7  = diff(F1,'N1') - W1;        % country 1's wage in terms of country 1's intermed. good
e8  = diff(F2,'N2') - W2;        % country 2's wage in terms of country 2's intermed. good
e9  = diff(F1,'K1') - R1;        % country 1's rental rate in terms of country 1's intermed. good
e10 = diff(F2,'K2') - R2;        % country 2's rental rate in terms of country 2's intermed. good
e11 = RX - (QA1/QA2);            % Law of one price
e12 = RX - (QB1/QB2);            % Law of one price


% RESCOURCE CONSTRAINTS AND PROD. FCT.:

e13 = ARM1 - (C1+X1+G1);                 % country 1 resource constraint, final good sector
e14 = ARM2 - (C2+X2+G2);                 % country 2 resource constraint, final good sector
e15 = A1+A2 - F1;                        % country 1 resource constraint, intermediate sector
e16 = B1+B2 - F2;                        % country 2 resource constraint, intermediate sector


% HOUSEHOLDS:

e17 =  diff(U1,'C1') - LAM1;                                                            % Marg. utility of cons., country 1
e18 =  diff(U2,'C2') - LAM2;                                                            % Marg. utility of cons., country 2      
e19 = -diff(U1,'N1') - (LAM1*QA1*W1);                                                   % Labor-leisure choice, country 1
e20 = -diff(U2,'N2') - (LAM2*QB2*W2);                                                   % Labor-leisure choice, country 2       
e21 = LAM1*(1+phiK*(K1p-K1)) - (betta*LAM1p*(1-delta+QA1p*R1p+phiK*(K1leadp-K1lead)));  % Capital EE, country 1
e22 = LAM2*(1+phiK*(K2p-K2)) - (betta*LAM2p*(1-delta+QB2p*R2p+phiK*(K2leadp-K2lead)));  % Capital EE, country 2
e23 = K1p - ((1-delta)*K1+X1);                                                          % Capital law of motion, country 1
e24 = K2p - ((1-delta)*K2+X2);                                                          % Capital law of motion, country 2
e25 =                       LAM1    *QQ - betta* LAM1p;                                 % Bond EE, country 1 (bond denominated in country 1's final good)
e26 = (1-quc*(log(Bp)/RX))*(LAM2/RX)*QQ - betta*(LAM2p/RXp);                            % Bond EE, country 2 (bond denominated in country 1's final good)

e33=C1+X1+G1+QQ*log(Bp) - (QA1*(R1*K1+W1*N1)+log(B)-phiK/2*(K1p-K1)^2);                 % Budget constraint, country 1
% % alternatively:
% e33=C2+X2+G2-QQ*log(Bp)/RX - (QB2*(R2*K2+W2*N2)-log(B)/RX-phiK/2*(K2p-K2)^2);           % Budget constraint, country 2

% EXOGENEOUS PROCESSES:

e27 = (Z1p) - (rho_Z1*Z1 +rho_Z2Z1*Z2 +(1-rho_Z1)*1    +(1-rho_Z2Z1)*1);
e28 = (Z2p) - (rho_Z2*Z2 +rho_Z1Z2*Z1 +(1-rho_Z2)*1    +(1-rho_Z1Z2)*1);
e34 = (G1p) - (rho_G1*G1 +rho_G2G1*G2 +(1-rho_G1)*G1bar+(1-rho_G2G1)*G2bar);
e35 = (G2p) - (rho_G2*G2 +rho_G1G2*G1 +(1-rho_G2)*G2bar+(1-rho_G1G2)*G1bar);


% AUXILIARY VARIABLES:

e36 = -K1lead + K1p;
e37 = -K2lead + K2p;

% OTHER VARIABLES OF INTEREST:

e29 = P        - (QB1/QA1);           % terms of trade, price of imports into country 1 relative to exports from country 1
e30 = GDP1     - (QA1*F1);            % country 1's GDP in units of country 1's final consumption good
e31 = IR1      - (B1/A1);             % country 1's import ratio at base year (const) prices
e32 = log(NX1) - (A2-P*B1)/GDP1;      % country 1's net Exports (as a fraction of GDP)

eta = [0         0       0        0;
       0         0       0        0;
       0         0       0        0;
       sigmaZ    0       0        0;
       0       sigmaZ    0        0; 
       0         0     sigmaG     0;
       0         0       0      sigmaG];
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE FUNCTION f, DEFINE CONTROLS AND STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22;e23;e24;e25;e26;e27;e28;e29;e30;e31;e32;e33;e34;e35;e36;e37];

% Define the vector of controls, y, and states, x
x  = [B  K1  K2  G1  G2  Z1  Z2];
xp = [Bp K1p K2p G1p G2p Z1p Z2p];
y  = [C1  C2  N1  N2  ARM1  ARM2  X1  X2  P  RX  A1  B1  A2  B2  W1  W2  R1  R2  LAM1  LAM2  QA1  QB1  QA2  QB2  QQ  IR1  GDP1  NX1  K1lead  K2lead];  
yp = [C1p C2p N1p N2p ARM1p ARM2p X1p X2p Pp RXp A1p B1p A2p B2p W1p W2p R1p R2p LAM1p LAM2p QA1p QB1p QA2p QB2p QQp IR1p GDP1p NX1p K1leadp K2leadp];  

%Make f a function of the logarithm of the state and control vector
f = subs(f, [x,y,xp,yp], exp([x,y,xp,yp]));

%Compute analytical derivatives of f
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp,approx);
