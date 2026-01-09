%%
% Sabău Alexandra-Denisa
%

clearvars
clc

%% Magic numbers
m = 6; 
n = 13; 

%% Process data (fixed, do not modify)
a1 = 2*(0.15+(m+n/20)/30)*(1000+n*300);
a2 = (1000+n*300);
b0 = (2.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [(-1)^n*(-m/10-rand(1)*m/5); (-1)^m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 20/a1/1e4; % fundamental step size
Tfin = 36/a1; % simulation duration

gain = 15;
umin = -gain; umax = gain; % input saturation
ymin = -b0*gain/1.8; ymax = b0*gain/1.8; % output saturation

whtn_pow_in = 1e-9*5*(((m-1)*5+n)/5)/2; % input white noise power and sampling time
whtn_Ts_in = Ts*3;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(9); % input quantizer (DAC)

whtn_pow_out = 1e-8*5*(((m-1)*8+n)/5)/2; % output white noise power and sampling time
whtn_Ts_out = Ts*5;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(9); % output quantizer (ADC)

u_op_region = -(m+n/5)/2; % operating point

%% Input setup (can be changed/replaced/deleted)
u0 = 0;     % fixed
ust = 4;  % must be modified (saturation) - pe asta o putem modifica 
t1 = 12/a1; % recommended 

%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_electricP1_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

plot(t,u,t,y)
shg

%% System identification
 i1 = 3651;
 i2 = 4886;
 i3 = 9890;
 i4 = 10980;
 
u0= mean(u(i1:i2));
ust= mean(u(i3:i4));
y0 = mean(y(i1:i2));
yst = mean(y(i3:i4));

K = (yst - y0)/(ust - u0)
%% partea reala a polilor
i5 = 5983;
i6 = 11910;

t_aux = t(i5:i6)
y_aux = abs(y(i5:i6) - yst);

figure 
plot(t_aux, y_aux )
i7 = 15;
i8 = 1239;
i9 = 2340;  % 3 maxime 

t_reg = t_aux([i7,i8,i9]);
y_reg = log(y_aux([i7,i8,i9]));

figure 
plot(t_reg, y_reg) % ar trebui sa fie o dreapta(pe cat se poate) 

A_reg = [sum(t_reg.^2), sum(t_reg) ; sum(t_reg), length(t_reg)];

b_reg = [sum(y_reg.*t_reg);
    sum(y_reg)];

theta = inv(A_reg)*b_reg
Re = theta(1)
%% partea imaginara  a polilor - prin diferenta dintre max si min, Im = wosc

i10 = 7240;
i11 = 8601;
Tosc =  2*(t(i11)-t(i10));
Im = 2*pi/ Tosc

%% zeta, wn 
wn = sqrt(Re^2 + Im^2)
zeta = -Re/wn

%% validarea

A = [0, 1 ; -wn^2, -2*zeta*wn];
B = [0; K*wn^2];
C = [1,0]
D = 0;
sys = ss(A,B,C,D);
ysim2 = lsim(sys,u, t,[y(1), 15]); %cond pt derivata , adica 15 
figure
plot(t,u,t,y,t,ysim2)
  %erori
J = 1/sqrt(length(t1))*norm(y-ysim2)
eMPN = norm(y-ysim2)/norm(y-mean(y))*100
  