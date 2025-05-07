% BKK95_param.m
% 
% March 2008, Katrin Rabitsch


function [omeg,sig,betta,quc,rho_Z1,rho_Z2,delta,gam,alfa,phiK,rho_Z1Z2,rho_Z2Z1,Z1bar,Z2bar,N1bar,N2bar,rho_G1,rho_G2,rho_G1G2,rho_G2G1,GY1bar,GY2bar,sigmaZ,sigmaG,eta]=BKK95_param;

betta=.99;          % discount factor
delta=.025;         % quarterly depreciation rate
gam=2;              % coefficient of relative risk aversion
sig=1.5;            % intratemporal elasticity of subsitution (between a and b goods)
omeg=.85;           % degree of home bias in final good production, implies import share of 15%
alfa=.36;           % Cobb-Douglas capital share
quc=.0005;          % portfolio adjustment cost (quadratic) parameter
phiK=0.5;           % capital adjustment cost (quadratic) parameter
% mue1=.34;         % parameter that pins down stst labor supply (such that N1=N1bar, etc.)
% mue2=.34;         % parameter that pins down stst labor supply (such that N2=N2bar, etc.)

rho_Z1=.906;        % persistence of technology shock, country 1
rho_Z2=rho_Z1;      % persistence of technology shock, country 2
rho_Z1Z2=.088;      % spillover of country 1 technology shock to country 2
rho_Z2Z1=rho_Z1Z2;  % spillover of country 1 technology shock to country 2

rho_G1=.906;        % persistence of government expenditure shock, country 1
rho_G2=rho_G1;      % persistence of government expenditure shock, country 2
rho_G1G2=.088;      % spillover of country 1 technology shock to country 2
rho_G2G1=rho_G1G2;  % spillover of country 1 technology shock to country 2

Z1bar=1;            % steady state technology level, country 1
Z2bar=1;            % steady state technology level, country 2
N1bar=.305;         % steady state labor supply, country 1 
N2bar=.305;         % steady state labor supply, country 2
GY1bar=0.2;         % steady state ratio of government expenditure to final output (Armington aggregator)
GY2bar=0.2;         % steady state ratio of government expenditure to final output (Armington aggregator)

sigmaZ=0.00852;
sigmaG=0.00852;


if sig==1; dummy_sig=1; else dummy_sig=0; end  % dummy, =1 if Armington aggregator becomes Cobb-Douglas

sigmaZ1Z2=0.00852;

eta=[zeros(3,4);eye(4)*sigmaZ1Z2];
