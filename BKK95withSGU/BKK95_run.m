% BKK95_run.m 
% 
% Two Country Two Good Model, following Backus, David K., Patrick J. Kehoe, and Finn E. Kydland [1995]. “International Business Cycles:
% Theory and Evidence,” in Thomas F. Cooley (ed.) Frontiers of Business Cycle Research, Princeton University Press, Princeton, 331-56.
%
% with slight modifications to original BKK95:
% -) incomplete markets (bond economy)
% -) capital adjustment costs (instead of time to build)
%
% Main file; calls the model files, gets solution, and plots impulse responses
%
% uses the following Matlab files:
% -) anal_deriv.m, gx_hx, gxx_hxx.m, gss_hss.m, num_eval.m, ir.m (by Stephanie Schmidt-Grohe and Martin Uribe)
% -) solab.m (by Paul Klein)
% -) csolve.m (by Christopher Sims)
%
% written by Katrin Rabitsch, March 2008
% 
%
% The linear solution is of the form
% y(t)   = gx x(t)
% x(t+1) = hx x(t) + e(t+1)
% 
% where x(t) is given by the vector of state variables   x = [B K1 K2 G1 G2 Z1 Z2]' and
% where y(t) is given by the vector of control variables y = [C1 C2 N1 N2 ARM1 ARM2 X1 X2 P RX A1 B1 A2 B2 W1 W2 R1 R2 LAM1 LAM2 QA1 QB1 QA2 QB2 QQ IR1 GDP1 NX1 K1lead K2lead]'
% 

global T dummy_sig shock


T      = 40;    % number of periods for plotting impulse responses
approx = 1;     % Order of approximation desired (==1: 1st order, ==2: 2nd order)
shock  = 5;     % ==1: government expenditure shock in country 1
                % ==2: government expenditure shock in country 2
                % ==3: productivity shock in country 1
                % ==4: productivity shock in country 2
                % ==5: plot impulse responses to all shocks

                
%Defines system of equations and takes derivatives
 [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f,eta] = BKK95_model(approx);
 if approx==1;
 anal_deriv_print2f('BKK95',fx,fxp,fy,fyp,f,eta);   
 elseif approx==2;
 anal_deriv_print2f('BKK95',fx,fxp,fy,fyp,f,eta, fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx);
 end
 
%Numerical evaluation: assign values to parameters and steady-state variables
 [B,K1,K2,Z1,Z2,C1,C2,N1,N2,ARM1,ARM2,X1,X2,P,RX,A1,B1,A2,B2,W1,W2,R1,R2,LAM1,LAM2,QA1,QB1,QA2,QB2,QQ,IR1,GDP1,NX1,Bp,K1p,K2p,Z1p,Z2p,C1p,C2p,N1p,N2p,ARM1p,ARM2p,X1p,X2p,Pp,RXp,A1p,B1p,A2p,B2p,W1p,W2p,R1p,R2p,LAM1p,LAM2p,QA1p,QB1p,QA2p,QB2p,QQp,IR1p,GDP1p,NX1p,G1,G2,G1p,G2p,K1lead,K2lead,K1leadp,K2leadp,mue1,mue2,G1bar,G2bar] = BKK95_stst;
 [omeg,sig,betta,quc,rho_Z1,rho_Z2,delta,gam,alfa,phiK,rho_Z1Z2,rho_Z2Z1,Z1bar,Z2bar,N1bar,N2bar,rho_G1,rho_G2,rho_G1G2,rho_G2G1,GY1bar,GY2bar,sigmaZ,sigmaG,eta] = BKK95_param;

%Obtain numerical derivatives of f
 BKK95_num_eval
   
%First-order approximation
 [gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp)
   
%Second-order approximation
 if approx==2;
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
 end;

% compute IR
 x0=zeros(size(hx,1),1);  
 if     shock==1 x0(end-3,1)=1; IR=ir(gx,hx,x0,T); BKK95_IR
 elseif shock==2 x0(end-2,1)=1; IR=ir(gx,hx,x0,T); BKK95_IR
 elseif shock==3 x0(end-1,1)=1; IR=ir(gx,hx,x0,T); BKK95_IR
 elseif shock==4 x0(end,1)=1;   IR=ir(gx,hx,x0,T); BKK95_IR
 end
 
 
 if shock==5
     for shock=1:4
         x0=zeros(size(hx,1),1); x0(end+shock-4,1)=1;  IR=ir(gx,hx,x0,T); BKK95_IR
     end
 end
     



