% BKK95_stst.m
% This file finds the non-stochastic steady state using csolve, calls file
% (BKK95_stst_equations.m) containing the system of equations
% 
% Katrin Rabitsch, March 2008

function [B,K1,K2,Z1,Z2,C1,C2,N1,N2,ARM1,ARM2,X1,X2,P,RX,A1,B1,A2,B2,W1,W2,R1,R2,LAM1,LAM2,QA1,QB1,QA2,QB2,QQ,IR1,GDP1,NX1,Bp,K1p,K2p,Z1p,Z2p,C1p,C2p,N1p,N2p,ARM1p,ARM2p,X1p,X2p,Pp,RXp,A1p,B1p,A2p,B2p,W1p,W2p,R1p,R2p,LAM1p,LAM2p,QA1p,QB1p,QA2p,QB2p,QQp,IR1p,GDP1p,NX1p,G1,G2,G1p,G2p,K1lead,K2lead,K1leadp,K2leadp,mue1,mue2,G1bar,G2bar]=BKK95_stst;

% assign parameters
[omeg,sig,betta,quc,rho_Z1,rho_Z2,delta,gam,alfa,phiK,rho_Z1Z2,rho_Z2Z1,Z1bar,Z2bar,N1bar,N2bar,rho_G1,rho_G2,rho_G1G2,rho_G2G1,GY1bar,GY2bar,sigmaZ,sigmaG,eta]=BKK95_param;

% vector of starting values
% x0=ones(27,1)*.5; 
x0=[0.6141;
    0.6141;
    0.2715;
    0.2715;
    1.1297;
    1.1297;
    0.2897;
    0.2897;
    0.0000;
    1.0000;
    0.9603;
    0.1695;
    0.1695;
    0.9603;
    2.3706;
    2.3706;
    0.0351;
    0.0351;
    0.6579;
    0.6579;
    1.0000;
    1.0000;
    1.0000;
    1.0000;
    0.9900;
   11.5867;
   11.5867];

% calling csolve
[SS,rc]=csolve(@BKK95_stst_equations,x0,[],0.00001,1000)

% displaying and storing output
fprintf('\n\n  STEADY STATE VALUES \n');
if rc==0;     fprintf('\n  (normal solution) \n\n') 
elseif rc==4; fprintf('\n  (WARNING: maximum number of iterations reached) \n\n')
else          fprintf('\n  (WARNING: no solution) \n\n') 
end

C1     = SS(1);           fprintf('     C1 = %8.5f\n',C1);      
C2     = SS(2);           fprintf('     C2 = %8.5f\n',C2);       
% N1     = SS(3);           fprintf('     N1 = %8.5f\n',N1);      
% N2     = SS(4);           fprintf('     N2 = %8.5f\n',N2);       
N1     =  N1bar;          fprintf('     N1 = %8.5f\n',N1);      
N2     =  N2bar;          fprintf('     N2 = %8.5f\n',N2);       
mue1   = SS(3);           fprintf('   mue1 = %8.5f\n',mue1);      
mue2   = SS(4);           fprintf('   mue2 = %8.5f\n',mue2);       
ARM1   = SS(5);           fprintf('   ARM1 = %8.5f\n',ARM1);      
ARM2   = SS(6);           fprintf('   ARM2 = %8.5f\n',ARM2);       
X1     = SS(7);           fprintf('     X1 = %8.5f\n',X1);      
X2     = SS(8);           fprintf('     X2 = %8.5f\n',X2);      
RX     = SS(10);          fprintf('     RX = %8.5f\n',RX);     
A1     = SS(11);          fprintf('     A1 = %8.5f\n',A1);     
B1     = SS(12);          fprintf('     B1 = %8.5f\n',B1);     
A2     = SS(13);          fprintf('     A2 = %8.5f\n',A2);               
B2     = SS(14);          fprintf('     B2 = %8.5f\n',B2);         
W1     = SS(15);          fprintf('     W1 = %8.5f\n',W1);        
W2     = SS(16);          fprintf('     W2 = %8.5f\n',W2);         
R1     = SS(17);          fprintf('     R1 = %8.5f\n',R1);        
R2     = SS(18);          fprintf('     R2 = %8.5f\n',R2);         
LAM1   = SS(19);          fprintf('   LAM1 = %8.5f\n',LAM1);   
LAM2   = SS(20);          fprintf('   LAM2 = %8.5f\n',LAM2);        
QA1    = SS(21);          fprintf('    QA1 = %8.5f\n',QA1);         
QB1    = SS(22);          fprintf('    QB1 = %8.5f\n',QB1);        
QA2    = SS(23);          fprintf('    QA2 = %8.5f\n',QA2);         
QB2    = SS(24);          fprintf('    QB2 = %8.5f\n',QB2);   
QQ     = SS(25);          fprintf('     QQ = %8.5f\n',QQ);         
B      = SS(9);           fprintf('      B = %8.5f\n',B);   
K1     = SS(26);          fprintf('     K1 = %8.5f\n',K1);        
K2     = SS(27);          fprintf('     K2 = %8.5f\n',K2);
K1lead = K1;
K2lead = K2;
P      = QB1/QA1;         fprintf('      P = %8.5f\n',P);      
Z1     = Z1bar;
Z2     = Z2bar;
G1bar  = GY1bar*ARM1;
G2bar  = GY2bar*ARM2;
F1     = Z1*K1^alfa*N1^(1-alfa);
F2     = Z2*K2^alfa*N2^(1-alfa);
G1     = G1bar;           fprintf('     G1 = %8.5f\n',G1);   
G2     = G2bar;           fprintf('     G2 = %8.5f\n',G2);         
IR1    = B1/A1;           fprintf('    IR1 = %8.5f\n',IR1);         
GDP1   = QA1*F1;          fprintf('   GDP1 = %8.5f\n',GDP1);   
NX1    = (A2-P*B1)/GDP1;  fprintf('    NX1 = %8.5f\n',NX1);         
G1     = G1bar;           fprintf('     G1 = %8.5f\n',G1);   
G2     = G2bar;           fprintf('     G2 = %8.5f\n',G2);         


% Applying logs and defining 'prime'-variables

C1   = log(C1);                     C1p    = C1;
C2   = log(C2);                     C2p    = C2;
N1   = log(N1);                     N1p    = N1;
N2   = log(N2);                     N2p    = N2;
ARM1 = log(ARM1);                   ARM1p  = ARM1;
ARM2 = log(ARM2);                   ARM2p  = ARM2;
X1   = log(X1);                     X1p    = X1;
X2   = log(X2);                     X2p    = X2;
P    = log(P);                      Pp     = P;
RX   = log(RX);                     RXp    = RX;
QQ   = log(QQ);                     QQp    = QQ;
A1   = log(A1);                     A1p    = A1;
A2   = log(A2);                     A2p    = A2;
B1   = log(B1);                     B1p    = B1;
B2   = log(B2);                     B2p    = B2;
W1   = log(W1);                     W1p    = W1;
W2   = log(W2);                     W2p    = W2;
R1   = log(R1);                     R1p    = R1;
R2   = log(R2);                     R2p    = R2;
LAM1 = log(LAM1);                   LAM1p  = LAM1;
LAM2 = log(LAM2);                   LAM2p  = LAM2;
QA1  = log(QA1);                    QA1p   = QA1;
QB1  = log(QB1);                    QB1p   = QB1;
QA2  = log(QA2);                    QA2p   = QA2;
QB2  = log(QB2);                    QB2p   = QB2;
IR1  = log(IR1);                    IR1p   = IR1;
GDP1 = log(GDP1);                   GDP1p  = GDP1;
Z1   = log(Z1);                     Z1p    = Z1;
Z2   = log(Z2);                     Z2p    = Z2;
G1   = log(G1);                     G1p    = G1;
G2   = log(G2);                     G2p    = G2;
K1   = log(K1);                     K1p    = K1;   
K2   = log(K2);                     K2p    = K2; 
K1lead = log(K1lead);               K1leadp = K1lead;   
K2lead = log(K2lead);               K2leadp = K2lead; 

% no logs since linear (not log-linear) approximation is taken for
% variables that are 0 at steady state and where percent deviations are not defined
NX1  = (NX1);                       NX1p  = NX1; 
B    = (B);                         Bp    = B;









%%
% BKK95_stst_equations.m
% defines the system of steady state equations
% 
% KR March 2008

function y=BKK95_stst_equations(x);

global dummy_sig

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assign parameters
[omeg,sig,betta,quc,rho_Z1,rho_Z2,delta,gam,alfa,phiK,rho_Z1Z2,rho_Z2Z1,Z1bar,Z2bar,N1bar,N2bar,rho_G1,rho_G2,rho_G1G2,rho_G2G1,GY1bar,GY2bar,sigmaZ,sigmaG,eta]=BKK95_param;

[rows,cols]=size(x);

j=1;
while j<=cols;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% variables in the system:

C1   = x(1,j);      C1p   = x(1,j);
C2   = x(2,j);      C2p   = x(2,j);
% N1   = x(3,j);     N1p    = x(3,j);
% N2   = x(4,j);     N2p    = x(4,j);
mue1 = x(3,j);   
mue2 = x(4,j);   
ARM1 = x(5,j);      ARM1p = x(5,j);
ARM2 = x(6,j);      ARM2p = x(6,j);
X1   = x(7,j);      
X2   = x(8,j);
B    = x(9,j);      Bp    = x(9,j);
RX   = x(10,j);     RXp   = x(10,j);
A1   = x(11,j);     A1p   = x(11,j);
B1   = x(12,j);     B1p   = x(12,j);
A2   = x(13,j);     A2p   = x(13,j);
B2   = x(14,j);     B2p   = x(14,j);
W1   = x(15,j);     W1p   = x(15,j);
W2   = x(16,j);     W2p   = x(16,j);
R1   = x(17,j);     R1p   = x(17,j);
R2   = x(18,j);     R2p   = x(18,j);
LAM1 = x(19,j);     LAM1p = x(19,j);
LAM2 = x(20,j);     LAM2p = x(20,j); 
QA1  = x(21,j);     QA1p  = x(21,j);   
QB1  = x(22,j);     QB1p  = x(22,j);   
QA2  = x(23,j);     QA2p  = x(23,j); 
QB2  = x(24,j);     QB2p  = x(24,j); 
QQ   = x(25,j);     QQp   = x(25,j);    
K1   = x(26,j);     K1p   = x(26,j);
K2   = x(27,j);     K2p   = x(27,j); 


Z1 = Z1bar; Z1p = Z1;
Z2 = Z2bar; Z2p = Z2;
N1 = N1bar; N1p = N1;
N2 = N2bar; N2p = N2;

F1  = Z1 *K1 ^alfa*N1 ^(1-alfa);
F1p = Z1p*K1p^alfa*N1p^(1-alfa);
F2  = Z2 *K2 ^alfa*N2 ^(1-alfa);
F2p = Z2p*K2p^alfa*N2p^(1-alfa);

% diff(U1,'C1')
% diff(U2,'C2')
% diff(U1,'N1')
% diff(U2,'N2')
% diff(F1,'N1')
% diff(F1,'K1')
% diff(F2,'N2')
% diff(F2,'K2')

U1_C1 =  (C1^mue1*(1-N1)^(1-mue1))^(1-gam)*mue1/C1;
U2_C2 =  (C2^mue2*(1-N2)^(1-mue2))^(1-gam)*mue2/C2;
U1_N1 = -(C1^mue1*(1-N1)^(1-mue1))^(1-gam)*(1-mue1)/(1-N1);
U2_N2 = -(C2^mue2*(1-N2)^(1-mue2))^(1-gam)*(1-mue2)/(1-N2);
F1_N1 = Z1*K1^alfa*N1^(1-alfa)*(1-alfa)/N1;
F2_N2 = Z2*K2^alfa*N2^(1-alfa)*(1-alfa)/N2;
F1_K1 = Z1*K1^alfa*alfa/K1*N1^(1-alfa);
F2_K2 = Z2*K2^alfa*alfa/K2*N2^(1-alfa);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SYSTEM OF STEADY STATE EQUATIONS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y(1,j) = A1-(   omeg *(QA1)^(-sig)*ARM1);  % country 1's demand for a good
y(2,j) = B2-(   omeg *(QB2)^(-sig)*ARM2);  % country 2's demand for a good
y(3,j) = B1-((1-omeg)*(QB1)^(-sig)*ARM1);  % country 1's demand for b good
y(4,j) = A2-((1-omeg)*(QA2)^(-sig)*ARM2);  % country 2's demand for b good

if dummy_sig==1;
    y(5,j) = ARM1 - (((A1^omeg)*(B1^(1-omeg)))/((omeg^omeg)*((1-omeg)^(1-omeg)))); % Armington aggregator, final output in country 1
    y(6,j) = ARM2 - (((A2^(1-omeg))*(B2^omeg))/((omeg^omeg)*((1-omeg)^(1-omeg)))); % Armington aggregator, final output in country 2  
else;
    y(5,j) = ARM1 - (((   omeg ^(1/sig))*(A1^((sig-1)/sig))+((1-omeg)^(1/sig))*(B1^((sig-1)/sig)))^(sig/(sig-1)));  % Armington aggregator, final output in country 1 
    y(6,j) = ARM2 - ((((1-omeg)^(1/sig))*(A2^((sig-1)/sig))+(   omeg ^(1/sig))*(B2^((sig-1)/sig)))^(sig/(sig-1)));  % Armington aggregator, final output in country 2     
end;

% INTERMEDIATE GOOD PRODUCERS:
y(7,j)  = F1_N1 - W1;                  % country 1's wage in terms of country 1's intermed. good
y(8,j)  = F1_K1 - R1;                  % country 1's rental rate in terms of country 1's intermed. good
y(9,j)  = F2_N2 - W2;                  % country 2's wage in terms of country 2's intermed. good
y(10,j) = F2_K2 - R2;                  % country 2's rental rate in terms of country 2's intermed. good
y(11,j) = RX - (QA1/QA2);              % Law of one price
y(12,j) = RX - (QB1/QB2);              % Law of one price

% RESCOURCE CONSTRAINTS AND PROD. FCT.:
y(13,j) = 1-(C1/ARM1+X1/ARM1+GY1bar);      % country 1 resource constraint, final good sector
y(14,j) = 1-(C2/ARM2+X2/ARM2+GY2bar);      % country 2 resource constraint, final good sector
y(15,j) = A1+A2-F1;                        % country 1 resource constraint, intermediate sector
y(16,j) = B1+B2-F2;                        % country 2 resource constraint, intermediate sector

% HOUSEHOLDS:
y(17,j) =  U1_C1 -  LAM1;                                   % Marg. utility of cons., country 1
y(18,j) =  U2_C2 -  LAM2;                                   % Marg. utility of cons., country 2
y(19,j) = -U1_N1 - (LAM1*QA1*W1);                           % Labor-leisure choice, country 1
y(20,j) = -U2_N2 - (LAM2*QB2*W2);                           % Labor-leisure choice, country 2
y(21,j) = 1 - (betta*(1-delta+R1p*QA1p));                   % Capital EE, country 1
y(22,j) = 1 - (betta*(1-delta+R2p*QB2p));                   % Capital EE, country 2
y(23,j) = K1p - ((1-delta)*K1+X1);                          % Capital law of motion, country 1
y(24,j) = K2p - ((1-delta)*K2+X2);                          % Capital law of motion, country 2
y(25,j) =                     QQ - (betta);                 % Bond EE, country 1 (bond denominated in country 1's intermed. good)
y(26,j) = (1-quc*1e3*(Bp/RX))*QQ - (betta);                 % Bond EE, country 2 (bond denominated in country 1's intermed. good)
y(27,j) = C1+X1+GY1bar*ARM1+QQ*Bp - (QA1*(R1*K1+W1*N1)+B);  % Budget constraint, country 1

j=j+1;
end;
