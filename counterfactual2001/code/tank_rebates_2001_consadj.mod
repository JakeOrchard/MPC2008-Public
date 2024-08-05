%--------------------------------------------------------------------------
% TANK_REBATES_2001_CONSADJ.MOD
%    PROGRAM TO EXPLORE 2001 REBATE COUNTERFACTUALS OF MICRO MPC IN NK DSGE
    
% By Valerie Ramey
%   Sept. 24, 2003 Version that raises g_ss and includes some as HH consumption
%      to deal with JPS ND scaling issue ("cgov_ss")
%   We want to target total cons = 67% of GDP and JPS ND = 53% of total cons in steady state
%   We don't want the rebate to be spent on non-JPS ND, so we take cgov_ss as an exog parameter

%   Revised Sept. 23, 2023 to add 6-month counterfactuals and correct calibration since JPS ND = 53% of cons
%   Revised March 29, 2023, 5:47 pm based on Johannes Wieland parameterization for 2008 program.
%   Comments out paramfile lines, sets phib=0.1 instead of 0.01.
%   However, still uses log utility (whereas Johannes' python programs use ies = 0.5)

% Conversion from quarterly to monthly calibration based on Sept. 20 and 26 2021 handwritten notes

% DETAILS OF PROGRAM

% Monthly calibration, uses path of 2001 rebates
% Medium-scale NK model
% Contains:
%   - rule-of-thumb consumers (TANK model)
%   - sticky prices and sticky wages with option of indexing of prices
%   - adjustment cost on investment rather than on capital
%   - variable utilization as in Schmidt-Grohe-Uribe
%   - possible shocks to transfers to ROT and optimizing households
%   - possible lags in paying off debt, allowing full deficit financing in short run
% Original program based on V. Ramey extensions of Sarah Zubairy program, Martin Uribe notes, 
%      Schmitt-Grohe-Uribe (sgu) (2005), Colciago (2011)

%--------------------------------------------------------------------------

%----------------------------------------------------------------
% 0. HOUSEKEEPING
%----------------------------------------------------------------

close all;
clc;

%----------------------------------------------------------------
% 1. VARIABLES AND PARAMETERS
%----------------------------------------------------------------

% VARIABLE DEFINITIONS

%   c = JPS consumption, aggregate
%   co = JPS consumption of optimizing agents
%   cr = JPS consumption of rule-of-thumb agents
%   cgov = part of govt spending that enters utility additively with JPS consumption

%   hd = labor hours demanded
%   h = hours, total supplied
%   ho = hours of optimizing agents
%   hr = hours of rule-of-thumb agents

%   iv = investment, total
%   ivo = investment of optimizing agents
%   k = capital stock, total, end of period
%   ko = capital of optimizing agents
%   q = Brainard-Tobin's q
%   u = capacity utilization

%   y = output
%   lambda = Lagrange multiplier on optimizing HH budget constraint
%   w = wage (real)
%   rk = rental price of capital (real)

%   g = government consumption purchases
%   b = government debt
%   t = total taxes
%   to = taxes on optimizing HH
%   tr = taxes on rule-of-thumb HH

%   s = price distortion due to price stickiness, s=1 in ss, >1 out of ss
%   sw = wage distortion due to wage stickiness, sw=1 in ss, >1 out of ss
%   r = gross nominal interest rate
%   pi = gross inflation rate, p/p(-1)
%   p = relative price of price adjusters to final goods price
%   mc = marginal cost
%   wstar = wage set by wage adjusters
%   x1 = auxiliary term for price stickiness
%   x2 = auxiliary term for price stickiness
%   f1 = auxiliary term for wage stickiness
%   f2 = auxiliary term for wage stickiness

%  ENDOGENOUS VARIABLES  (34 variables)

var co cr c;          % JPS consumption, aggregate JPS consumption
var ho hr h hd;       % hours 
var ko k ivo iv q u;  % investment and capital stock, utilization
var y;                % output
var lambda w rk;      % Lagrange multiplier on HH budget constraint, wages, rental rate of capital
var p mc muw;         % relative price, marginal cost, wage markup 
var s sw wstar ;      % price distortion( due to price stickiness), wage distortion, reset wage 
var x1 x2 f1 f2;      % auxiliary terms for price and wage stickiness
var r pi;             % gross nominal interest rate,inflation
var t tr to b;        % total taxes, taxes on ROT HH, taxes on optimizing HH, govt debt  

% EXOGENOUS VARIABLES (1 variable)

varexo e;        % shock to fiscal
 

%----------------------------------------------------------------
% 2. PARAMETERS AND CALIBRATION
%----------------------------------------------------------------

% PARAMETERS

parameters beta;                     % discount factor
parameters delta;                    % depreciation rate on private capital
parameters kappa;                    % investment adjustment cost
parameters delta1;                   % parameter of utilization cost
parameters delta2;                   % parameter of utilization cost
parameters ubar;                     % steady-state capital utilization
parameters alpha;                    % capital share
parameters gamma;                    % rule-of-thumb share
parameters theta;                    % price stickiness
parameters thetaw;                   % wage stickiness
parameters chi;                      % degree of wage indexing
parameters phipi phigap;             % Taylor rule parameters
parameters phib;                     % fiscal financing parameter
parameters eps;                      % elasticity of substitution of goods varieties
parameters epsw;                     % elasticity of substitution of labor varieties
parameters muw_ss;                   % steady-state wage markup
parameters mup_ss;                   % steady-state price markup
parameters phi;                      % utility parameter (1/phi is Frisch elasticity)
parameters nu;                       % weight on disutility of labor
parameters rhor;                     % Inertia in the Taylor Rule

% Steady state parameters

parameters co_ss cr_ss c_ss cy_ss cgovy_ss cgov_ss;
parameters ho_ss hr_ss h_ss hd_ss ko_ss k_ss ivo_ss iv_ss q_ss u_ss;
parameters lambda_ss w_ss rk_ss s_ss sw_ss r_ss pi_ss y_ss;
parameters t_ss tr_ss to_ss p_ss b_ss mc_ss wstar_ss x1_ss x2_ss f1_ss f2_ss g_ss gy_ss;


delta    = 1-(1-.015)^(1/3); % monthly depreciation rate of private capital (annual = 6%, quarterly 1.5%)
kappa    = 40;               % investment adj. cost, monthly value that approx IRF of 5.2 for quarterly from Leeper et al.
phi      = 1;                % 1/phi is Frisch elasticity) GLV =  0.2
alpha    = 0.36;             % capital share; GLV = 1/3
mup_ss   = 1.2;              % steady state price markup
muw_ss   = 1.2;              % wage markup
theta    = 0.9167;           % price stickiness, avg. duration of a price = 1/(1-theta). monthly equiv to theta = 0.75 at quarterly, 0 implies no stickiness; GLV = 0.75
thetaw   = 0.9167;           % wage stickiness, monthly equiv to theta = 0.75 at quarterly, 0 implies no stickiness
chi      = 0;                % wage indexing, , 0=no, 1=full
gamma    = 0.375;            % rule-of-thumb share: 3-month = 0.375 from JPS Table 2C, 6-month=0.66 from Table 4A (2SLS) (GLV = 0.5) 
beta     = 0.99^(1/3);       % monthly discount factor
phib     = 0.1;              % fiscal rule parameters, our quarterly = 0.1, GLV quarterly = 0.33 
phipi    = 1.5;              % Taylor rule parameter on inflation
phigap   = 1/12;             % Taylor rule parameter on output gap, 1 for annual in balanced-approach rule, divide by 12 for monthly
rhor     = 0.85^(1/3);       % Inertia in Taylor rule, 0 = no inertia, quarterly rhor = monthly rhor^3, .85 from Fed website
eps      = mup_ss/(mup_ss-1);% elasticity of subsitution of goods varieties
epsw     = muw_ss/(muw_ss-1);% elasticity of substitution of labor varieties
pi_ss    = 1;                % steady state inflation
nu       = 1;                % weight on disutility of labor
cy_ss    = 0.35;             % JPS ND consumption fraction of GDP is 35% in 2001 ****
delta1   = 1/beta-1+delta;   % utiliz. cost linear term, so u = ubar in SS
delta2   = 2*delta1;         % utiliz. cost quadratic term, SGU 2.02 based on ACEL, Leeper et al. approx. 2
ubar     = 1;                % steady-state utilization

%----------------------------------------------------------------
% 3. STEADY-STATE VALUES OF ENDOGENOUS VARIABLES
%----------------------------------------------------------------

% Note that the order of these statements matters.

s_ss = 1;
p_ss = 1;
mc_ss = 1/mup_ss;
rk_ss = 1/beta -1+delta;
sw_ss = 1;
u_ss = ubar;

ky_ss = alpha/((1/beta - 1 + delta)*mup_ss); 
iy_ss = delta*ky_ss;
cgovy_ss = 0.67 - cy_ss;  % Set the cons part of gov so that total consumption is 67% of GDP
gy_ss = 1 - iy_ss - cy_ss; % Total gov spending, including cgov
h_ss = ((1-alpha)/(mup_ss*muw_ss*nu*(cy_ss + cgovy_ss)))^(1/(1+phi)); 
hd_ss = h_ss;
y_ss = ((ky_ss^alpha)*(hd_ss^(1-alpha)))^(1/(1-alpha)); 
k_ss = y_ss* ky_ss;
w_ss = (k_ss/hd_ss)*rk_ss*(1-alpha)/alpha;
wstar_ss = w_ss;
c_ss = cy_ss*y_ss; 
cgov_ss = cgovy_ss*y_ss;
ho_ss = h_ss;
hr_ss = h_ss;
ko_ss = k_ss/(1-gamma);
iv_ss = delta*k_ss;
ivo_ss = iv_ss/(1-gamma);
q_ss = 1;
co_ss = c_ss;
cr_ss = c_ss;
g_ss = y_ss-iv_ss-c_ss;
tr_ss = w_ss*hr_ss-cr_ss;  % ***
r_ss = 1/beta;
b_ss = 0;
to_ss = (g_ss-gamma*tr_ss)/(1-gamma);
t_ss = (1-gamma)*to_ss + gamma*tr_ss;
x1_ss = y_ss/(1-theta*beta);
x2_ss = x1_ss;
lambda_ss = 1/(co_ss+cgov_ss);
f1_ss = w_ss*lambda_ss*hd_ss*muw_ss^(-1)/(1-thetaw*beta);
f2_ss = f1_ss;

%----------------------------------------------------------------
% MODEL EQUATIONS   (35 equations)
%----------------------------------------------------------------

% UTILITY FUNCTION FOR REFERENCE
%util = gamma*(ln(cr + cgov_ss) - nu*(hr^(1+phi))/(1+phi)) + (1-gamma)*(ln(co + cgov_ss) - nu*(ho^(1+phi))/(1+phi));

model; 

%----------------------------------------------------------------
% Expressions for use in model equations
%----------------------------------------------------------------

#adj = (kappa/2)*(ivo/ivo(-1) - 1)^2;
#adjdif = kappa*(ivo/ivo(-1) - 1);
#adjdifp = kappa*(ivo(+1)/ivo - 1);
#au = delta + delta1*(u - ubar) + (delta2/2)*(u - ubar)^2;
#aup = delta + delta1*(u(+1 )- ubar) + (delta2/2)*(u(+1) - ubar)^2;
#audif = delta1 + delta2*(u-ubar);


%----------------------------------------------------------------
% Model equations
%----------------------------------------------------------------

% 1) Capital equation of motion

  ko = (1-au)*ko(-1) + ivo*(1-adj);

% 2) SHADOW VALUE OF WEALTH

  lambda = 1/(co + cgov_ss);

% 3) CONSUMPTION EULER EQUATION WITH TOBIN'S Q 

 lambda*q = beta*lambda(+1)*(rk(+1)*u(+1) + q(+1)*(1-aup));

% 4) INVESTMENT ADJUSTMENT COST

  lambda = lambda*q*(1 - adj - (ivo/ivo(-1))*adjdif) + beta*lambda(+1)*q(+1)*(ivo(+1)/ivo)^2*adjdifp;

% 5) OPTIMAL CAPITAL UTILIZATION

  rk = q*audif;

% 6) RULE-OF-THUMB CONSUMPTION - Choose contemporaneous or with lags ***

%   cr = w*hr - tr; // contemporaneous

%  Use this for 3-month moving average case and gamma = 0.375

 cr = (w*hr - tr + w(-1)*hr(-1) - tr(-1) + w(-2)*hr(-2) - tr(-2))/3;

% Use this for 6-month moving average case and gamma = 0.66

%   cr = (0.387/0.66)*(w*hr - tr + w(-1)*hr(-1) - tr(-1) + w(-2)*hr(-2) - tr(-2))/3
%        + (1 - 0.387/0.66)*(w(-3)*hr(-3) - tr(-3) + w(-4)*hr(-4) - tr(-4) + w(-5)*hr(-5) - tr(-5))/3; 
 
% 7) MRS CONDITION

  muw = w*(gamma/((cr+cgov_ss)*nu*hr^phi) + (1-gamma)/((co+cgov_ss)*nu*ho^phi));

% 8) EQUALITY OF HOURS ACROSS TYPES

  ho = hr;

% 9)-11) SGU RECURSIVE EQS FOR OPTIMAL WAGE SETTING - TANK VERSION (substitutes gamma/cr + (1-gamma)/co for lambda in f1)

%  f1 = (epsw-1)/epsw*wstar*(gamma/cr + (1-gamma)/co)*(w/wstar)^epsw*hd + thetaw*beta*(pi(+1)/pi^chi)^(epsw-1)*(wstar(+1)/wstar)^(epsw-1)*f1(+1);  
f1 = (epsw-1)/epsw*wstar*(gamma/(cr+cgov_ss) + (1-gamma)/(co+cgov_ss))*(w/wstar)^epsw*hd + thetaw*beta*(pi(+1)/pi^chi)^(epsw-1)*(wstar(+1)/wstar)^(epsw-1)*f1(+1);  
  f2 = nu*h^phi*(w/wstar)^epsw*hd + thetaw*beta*(pi(+1)/pi^chi)^epsw*(wstar(+1)/wstar)^epsw*f2(+1); 
  f1 = f2;

% 12) - 13) DISTORTIONS DUE TO STICKY WAGES

  h = sw*hd;
  sw = (1-thetaw)*(wstar/w)^(-epsw) + thetaw*(w(-1)/w)^(-epsw)*(pi/pi(-1)^chi)^epsw*sw(-1);

% 14) AGGREGATE WAGES

  w^(1-epsw) = (1-thetaw)*wstar^(1-epsw) + thetaw*w(-1)^(1-epsw)*(pi(-1)^chi/pi)^(1-epsw);

% 15) FOC FOR FIRM OPTIMAL CAPITAL-LABOR RATIO

  u*k(-1)/hd = alpha/(1-alpha)*w/rk; 

% 16) REAL MARGINAL COST

  mc = rk^alpha*w^(1-alpha)*alpha^(-alpha)*(1-alpha)^(-1+alpha); 

% 17) - 19) URIBE RECURSIVE EQS FOR OPTIMAL PRICE SETTING (SGU 2005 SEEM DIFFERENT)

  x1 = mup_ss*p^(-eps)*y*mc + theta*beta*lambda(+1)/lambda*x1(+1);
  x2 = p^(1-eps)*y + theta*beta*lambda(+1)/lambda*p/p(+1)*x2(+1)/pi(+1);
  x1 = x2;

% 20) AGGREGATE INFLATION

  1 = (1-theta)*p^(1-eps) + theta*pi^(eps-1);

% 21) AGGREGATE JPS CONSUMPTION

  c = gamma*cr + (1-gamma)*co;

% 22) AGGREGATE HOURS

  h = gamma*hr + (1-gamma)*ho;

% 23) AGGREGATE INVESTMENT

  ivo=iv/(1-gamma);

% 24) AGGREGATE CAPITAL

  ko=k/(1-gamma);

% 25) AGGREGATE TAXES

  t = gamma*tr + (1-gamma)*to;

% 26) RESOURCE CONSTRAINT

  y = c + iv + g_ss;

% 27) PRODUCTION FUNCTION (S>=1 CAPTURES DISTORTION DUE TO STICKY PRICES)

  s*y = (u*k(-1))^alpha*hd^(1-alpha); 

% 28) EQUATION OF MOTION FOR STICKY PRICE DISTORTION

  s = (1-theta)*p^(-eps) + theta*pi^eps*s(-1);

% 29) GROSS NOMINAL INTEREST RATE

  lambda = r*beta*lambda(+1)/pi(+1);

% 30) TAYLOR RULE 

  r = rhor*r(-1) + (1-rhor)*(r_ss + phipi*(pi-pi_ss) + phigap*(y/y_ss - 1)); % ****

%  (r-r_ss) = phipi*(pi-pi_ss) + phigap*(y/y_ss - 1);

% 31) GOVERNMENT BUDGET CONSTRAINT

  g_ss = b/r - b(-1)/pi + t;

% 32) - 33) TAX RULES - ASSUMES THEY DON'T BEGIN RAISING LUMP-SUM TAXES TO PAY OFF THE DEFICIT FOR MANY MONTHS

  t = t_ss + phib*(b(-12) - b_ss) - e; // ***
  to = to_ss + tr - tr_ss;  // GLV implicit assumption (according to Uribe notes)

end;

%----------------------------------------------------------------
% STEADY-STATE
%----------------------------------------------------------------

initval;
  e = 0;
end;

steady_state_model;
  co = co_ss;
  cr = cr_ss;
  c = c_ss;
  ho = ho_ss;
  hr = hr_ss;
  h = h_ss;
  k = k_ss;
  ko = k_ss/(1-gamma);
  ivo = iv_ss/(1-gamma);
  iv = iv_ss;
  q = q_ss;
  u = u_ss;
  lambda = lambda_ss;
  w = w_ss;
  rk = rk_ss;
  s = s_ss;
  r = r_ss;
  pi = pi_ss;
  y = y_ss;
  to = (g_ss-gamma*tr_ss)/(1-gamma);
  tr = tr_ss;
  t = (1-gamma)*to_ss + gamma*tr_ss;
  p = p_ss;
  b = b_ss;
  mc = mc_ss; 
  x1 = x1_ss; 
  x2 = x2_ss; 
  wstar = wstar_ss;
  f1 = f1_ss;
  f2 = f2_ss;
  muw = muw_ss;
  hd = hd_ss;
  sw = sw_ss;
end;

 check;
 steady; 

% The following is calibrated to actual path of rebates in 2001 (2001q1 is y_ss)
shocks;
 var e;
 periods 2 3 4 5 6 7;
 values (0.008*y_ss) (0.020*y_ss) (0.013*y_ss) (0.0002*y_ss) (0.0005*y_ss) (0.0002*y_ss);  // fraction of GDP
% values (0.022*c_ss) (0.056*c_ss) (0.035*c_ss) (0.0007*c_ss) (0.001*c_ss) (0.0007*c_ss);  // fraction of JPS ND
end;

perfect_foresight_setup(periods=1200);
perfect_foresight_solver;

mpc = -sum(c(1:13) - c(1)) / sum(tr(1:13) - tr(1));
disp(sprintf('GE MPC is %.3f.',mpc))

% Copy the following A matrix to baseline_2001.irfs.xlsx

A = [oo_.exo_simul, y, c, co, cr, t, to, tr, h, ho, hr, w, iv, u, rk, r, pi, k, b] ;
% A = [oo_.exo_simul, oo_.endo_simul'];
% A = [oo_.exo_simul, g, y, c, co, cr, t, to, tr, h, ho, hr, w, iv, u, rk, r, pi, k, b] ;
%A = [oo_.exo_simul, g, y, c, co, cr, t, to, tr, h, ho, hr, w] ;


