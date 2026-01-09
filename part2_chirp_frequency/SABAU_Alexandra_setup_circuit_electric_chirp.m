%%
% Nume si prenume: Sabău Alexandra
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 6; %4
n = 13;  %32  

%% Process data (fixed, do not modify)
a1 = 2*(0.15+(m+n/20)/30)*(1000+n*300);
a2 = (1000+n*300);
b0 = (2.2+m+n)/5.5;

rng(m+10*n)
x0_slx = [(-1)^n*(-m/10-rand(1)*m/5); (-1)^m*(n/20+rand(1)*n/100)];

%% Experiment setup (fixed, do not modify)
Ts = 20/a1/1e4; % fundamental step size
Tfin = 36/a1*10; % simulation duration

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
wf = 4646 %modific; % ~wosc (foarte apropiat de wn care ar trebui de fapt, usor de citit din step 
fmin = wf/2/pi/10;
fmax = wf/2/pi*4;
Ain = 1.25
%% Data acquisition (use t, u, y to perform system identification)
out = sim("circuit_electric_R2022b_chirp.slx");

t = out.tout;
u = out.u;
y = out.y;


%filtru
 y=sgolayfilt(y,1,41);
 u=sgolayfilt(u,1,41);

plot(t,u,t,y)
shg

%% System identification
%maxim si min de la rosu (de unde incepe sinusul)

Ay=(-11.39+21.015)/2; 
Au = 1.25;
K = Ay/Au

%%

%pt punctul de rezonanta -max de la rosu 
% valorile de la y (ala de sus - ala de jos )
Ayr = (-15.3475+17.7407)/2
Aur = 1.25;
Mr = Ayr/Aur

%pulsatia de la rezonanta 
wr = pi/(0.0216718-0.0209344) % valorile de la x 
r = roots([4*Mr^2, 0, -4*Mr^2, 0,K^2])  % luam rad poz si pe cea mai mica . trebuie sa fie mai mica decat sqrt(2)/2
zeta = 0.3781;
wn =wr/sqrt(1-2*zeta^2);


%% validare 
 % cu k incepem ca sa vedem daca e centrat 
%K =3.85
%zeta = 0.3850
A = [0, 1 ; -wn^2, -2*zeta*wn];
B = [0; K*wn^2];
C = [1,0];
D = 0;
sys = ss(A,B,C,D);
ysim2 = lsim(sys,u, t,[y(1), +19000]); %cond pt derivata
figure
plot(t,u,t,y,t,ysim2)
  %erori
J = 1/sqrt(length(t))*norm(y-ysim2)
 eMPN = norm(y-ysim2)/norm(y-mean(y))*100


