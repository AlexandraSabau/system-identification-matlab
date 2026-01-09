%%
% Nume si prenume: Sabău Alexandra-Denisa
%

clearvars
clc

%% Magic numbers (replace with received numbers)
m = 6;
n = 13;

%% Process data and experiment setup (fixed, do not modify)
u_star = 0.15+n*0.045; % trapeze + PRBS amplitudes
delta = 0.02;
delta_spab = 0.015;

E = 12;  % converter source voltage

umin = 0; umax = 0.98; % input saturation
assert(u_star < umax-0.1)
ymin = 0; ymax = 1/(1-u_star)*E*2; % output saturation

% circuit components + parasitic terms
R = 15;
rL = 10e-3;
rC = 0.2;
rDS1 = 0.01;
rDS2 = 0.01;
Cv = 600e-6/3*m;
Lv = 40e-3*3/m;

% (iL0,uC0)
rng(m+10*n)
x0_slx = [(-1)^(n+1)*E/R,E/3/(1-u_star)];

Ts = 1e-5*(1+2*(u_star-0.15)/u_star); % fundamental step size
Ts = round(Ts*1e6)/1e6;

% input white noise power and sampling time
whtn_pow_in = 1e-11*(Ts*1e4)/2; 
whtn_Ts_in = Ts*2;
whtn_seed_in = 23341+m+2*n;
q_in = (umax-umin)/pow2(11); % input quantizer (DAC)

% output white noise power and sampling time
whtn_pow_out = 1e-7*E*(Ts*1e4/50)*(1+(50*u_star)*(u_star-0.15))/3; 
whtn_Ts_out = Ts*2;
whtn_seed_out = 23342-m-2*n;
q_out = (ymax-ymin)/pow2(11); % output quantizer (ADC)

meas_rep = 13+ceil(n/2); % data acquisition hardware sampling limitation

%% Input setup (can be changed/replaced/deleted)
%t1 = 0.3;
t1 = 0.2;
tr = 0.05 *5/2;
%tr = 0.02 * 3 % daca inmultesti cu 10 tr o sa fie iesirea deja va fi raspuns la treapta, prea mult , daca /10 va nu vei vedea cat vrei sa vezi tr
% inmultim ca sa prindem putin dintr o perioada % timp de urvcare * 2 sau 3 datorita oscilatiilor puternice 
N = 5; % la alegere - > 2^4 - 1 
% prima parte ramane fixa, primul bloc, dar dupa voi avea variatiuni de
% blocuri mai inguste 
p = round(tr/N/Ts);   % p = prescaler 

DeltaT = p*(2^N)*Ts*2;  
% poti sa i dai mai mult decat perioada asta , dar nu mai putin 
[input_LUT_dSPACE,Tfin] = generate_input_signal(Ts,t1,DeltaT, N,p,u_star,delta,delta_spab);

%LUT look up table - lista de perechi x, y => u, y 
% trebuie 2^4 - 1 esantioane 


%% Data acquisition (use t, u, y to perform system identification)
out = sim("convertor_R2022b.slx");

t = out.tout;
u = out.u;
y = out.y;

subplot(211)
plot(t,u)
subplot(212)
plot(t,y)

%% System identification
% nu dau mare factorul de umplere 
% cat de mult sa dureze un interval =p 
%0.3, 0.5 t1 pana se stabilizeaza iesirea 
%curs 5 , calibrare SPAB


% de sus pe margini 


i1 = 42702;
i2 = 113809;

%partea de jos de pe margini 


i3 = 146015;
i4 = 216399;

Nr = 19;

t_id = t(i1:Nr:i4);
u_id = u(i1:Nr:i4);
u_id = u_id-mean(u_id);
% u_id = detrend(u_id)
y_id = y(i1:Nr:i4);
 y_id = y_id-mean(y_id);
% y_id = detrend(y_id)

t_vd = t(i3:Nr:i4);
u_vd = u(i3:Nr:i4);
 u_vd = u_vd-mean(u_vd);
% u_vd = detrend(u_vd)

y_vd = y(i3:Nr:i4);
 y_vd = y_vd-mean(y_vd);
% y_vd = detrend(y_vd)

%afisare subsemnale 
figure

subplot(221)
plot(t_id, u_id)

subplot(223)
plot(t_id, y_id)

subplot(222)
plot(t_vd, u_vd)

subplot(224)
plot(t_vd, y_vd)

dat_id = iddata(y_id, u_id,t_id(2)-t_id(1));
dat_vd = iddata(y_vd,u_vd,t_vd(2)-t_vd(1));

%% MODEL ARX
model_arx = arx(dat_id,[2,2,1])
figure, resid(model_arx,dat_vd) %e valid pt ca primele 4 esant sunt in banda, 4 pt ca atatia coef am 
figure, compare(model_arx,dat_vd)
zpk(model_arx)

%% model armax - e bun -trece autocorelație
%curs 7 
model_arxmax= armax(dat_id,[5,5,5,1]) 
figure, resid(model_arxmax,dat_vd) 
figure, compare(model_arxmax,dat_vd)
zpk(model_arxmax)
%% OE
%validat pt xcorr - trece interc, dar fit foarte mic 
model_oe = oe(dat_id,[8,1,2])
figure, resid(model_oe,dat_vd) 
figure, compare(model_oe,dat_vd)

%% BJ
model_bj = bj(dat_id,[2 2 2 8 2])
figure, resid(model_bj, dat_vd)
figure, compare(model_bj, dat_vd)
zpk(model_bj)

%% iv4- a trecut 
model_iv4 = iv4(dat_id,[2 2 1])
figure, resid(model_iv4, dat_vd)
figure, compare(model_iv4, dat_vd)
pole(model_iv4)
dcgain(model_iv4)

%% cu spatiul starilor 
model_n4sid = n4sid(dat_id,14)
%model_n4sid = n4sid(dat_id,2)
figure, resid(model_n4sid,dat_vd) 
figure, compare(model_n4sid,dat_vd)
zpk(model_n4sid)
%%
model_ssest = ssest(dat_id,3) 
%model_n4sid = n4sid(dat_id,2)
figure, resid(model_ssest,dat_vd) 
figure, compare(model_ssest,dat_vd)
zpk(model_ssest) % 1 zero , 2 poli complex conjug, cu - pt ca avem faza neminima 
model_ssest.Ts
