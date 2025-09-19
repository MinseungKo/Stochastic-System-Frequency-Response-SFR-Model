clear all
close all
clc

%% Import Data
rng('shuffle')

PFile='loadbus_3_500.xlsx'; %% PSSE Output file
vsens = csvread('volt_sensitivity.csv'); %% Voltage sensitivity obtained by reactive power perturbation
H = csvread('voltH.csv');

data = xlsread(PFile);
data = data(5:end, :);

time = data(:, 1);
Ps = data(:, 2:11);
Freq = data(:, 12:50);
Volt = data(:, 51:end);

dt = time(2,1);
T = time(end,1);
N = round(T/dt)+1; % Total solve steps
sims = 1000; % Total simulation runs

%% Parameter Settings
cont_size = 500; % Contingency size (Base case: 500 MW)

% Governor and turbine model parameters
Rg = 0.05;            % Governor droop
Tg = 0.5;             % Governor time constant
Trh = 7;              % Reheat time constant
Fh = 0.3;             % Fraction of reheat power

% GFL Coefficients
RGFL = 0.045;
Tconv = 0.02; % Combined with measurement & activation delays
KGFL = 10;

% GFM Coefficients
Ome_b = 2*pi*60;
K_GFM = 3;
Droop_GFM = 0.045;
H_gfm = 4;
D_gfm = 2;

% FFR Coefficients
Rffr = 0.05;
Tffr_dev = 0.02;
Tffr_roc = 0.02;

%% Generator & System Informations
% In our case, Bus 31~39 are generator buses.
% Bus 31, 32, 33, 36, 38, 39: Synchronous generators
% Bus 34, 35, 37: Large IBRs

Basepow = 100;
Sbase = sum(Ps(1,:))*Basepow; %MW
SG_num = [2, 3, 4, 7, 9, 10];
Genbase = sum(Ps(1,SG_num))*Basepow;
GFLbase = Sbase-Genbase;

Hs = [2.3, 4.9, 4.3, 5.2, 3.15, 5.0]; % 31, 32, 33, 36, 38, 39
freqsys = sum(Hs.*Freq(:, SG_num),2)/sum(Hs);

Hsys = 0;
for i=1:6
    Hsys = Hsys+Hs(i)*Ps(1,SG_num(i))*100;
end
Dsys=0;
Hsys = Hsys/Sbase;

%% Dynamic Equations: System, SG, and Contingency
% System Dynamics
Num = tf(1, [2*Hsys, Dsys]);
[Anum, Bnum, Cnum, Dnum] = ssdata(Num);

% Generator
Gen = tf(1, [Tg, 1]) * tf([Fh * Trh, 1], [Trh, 1]) * (1 / Rg);
% AGC = tf(5, [1, 0]); 
% Gen1 = parallel(Gen, AGC); % Gen1 can be used instead for SG with AGC
[Agen, Bgen, Cgen, Dgen] = ssdata(Gen);

% Load Contingency
P_m = @(t) 0 * (t < 1) + (-cont_size/Sbase) * (t >= 1); % P_a is 0 until t = 1, then -0.1

%% Data for Voltage Estimation
% Initial load params
ini_load_v = [1.0345, 1.0221, 1.0102, 1.0098, 1.0109, 1.024, 1.0374, ...
                1.0355, 0.9941, 1.0359, 1.0474,1.0434,1.0586,1.0547,1.0428,1.0515,1.0518,0.982,1.03];
ini_load_p = [322, 400, 233.8, 522, 7.5, 220, 329.4, 158, 528, 274, 274.5, 208.6, 224, 139, 181, ...
                206,183.5, 9.2, 1104];
allbusnum = [30, 31, 32, 33, 34, 35, 36, 37, 38, 39];
load_bus = [3,4,7,8,12,15,16,18,20,21,23,24,25,26,27,28,29];
gen_load_bus = [2, 10];

for i=1:size(Volt,2)
    dvmin(i) = min(Volt(:,i)-Volt(1,i));
    [vmin(i), minind] = min(Volt(:,i)./Volt(1,i)-1);   
    minfreq = freqsys(minind);
end
vfrat = vmin/minfreq/60;

genvolts = Volt(:, allbusnum);
dgenvolts = genvolts-genvolts(1,:);
ratios = zeros(size(dgenvolts,2), size(genvolts,1));

% Conventional Method
for i=1:size(genvolts,1)
    conventional(i,:) = vsens(:, 1:29)'*dgenvolts(i,:)';
end

% Proposed Method
H_L = H(1:29, :);
H_G = H(30:end, :);

[u, s, v] = svd(H);

% est_Hl = u_l*s*v';
% est_Hg = u_g*s*v';
% est_H = u*s*v';

u_l = u(1:29, 1:10);
u_g = u(30:end, 1:10);
s = s(1:10, :);

for i=1:size(genvolts,1)
    dV = dgenvolts(i,:);
    dQest(i,:) = pinv(H_G)*dV';
    hat_VG = H_G*dQest(i,:)';
    residual(i,:) = dV-hat_VG';

    x = v'*dQest(i,:)';
    y = u_g'*dV';

    % Method1
    e_t = y - s*x;
    ds_t = e_t./x;
    s2 = s+diag(ds_t);

    H_Gt = u_g*s2*v';
    H_Lt = u_l*s2*v';
    S_LGt = H_Lt*pinv(H_Gt);
    proposed(i,:) = S_LGt*dV';
end
proposed(1,:)=0;
loadvolts = conventional+Volt(1, 1:29);
loadcorrvolts = proposed+Volt(1, 1:29);

load_bus = [3,4,7,8,12,15,16,18,20,21,23,24,25,26,27,28,29];
gen_load_bus = [2, 10];

figure(20)
minvoltcomp = [min(Volt(1:361,1:29)-Volt(1,1:29),[],1)', min(loadvolts(1:361,1:29)-loadvolts(1,1:29),[],1)', min(loadcorrvolts(1:361,1:29)-loadcorrvolts(1,1:29),[],1)'];
bar(minvoltcomp)

err_prop = zeros(size(load_bus));
err_conv = zeros(size(load_bus));
figure(21)
for i=1:size(load_bus,2)
    a = load_bus(i);
    plot(time, Volt(:,a), 'LineWidth', 2, 'Color','k')
    hold on
    plot(time, loadvolts(:,a), 'b--', 'LineWidth', 2)
    plot(time, loadcorrvolts(:,a), 'r:', 'LineWidth', 2)

    err_prop(1,i) = sum(abs(Volt(:,a)-loadcorrvolts(:,a)));
    err_conv(1,i) = sum(abs(Volt(:,a)-loadvolts(:,a)));
end

DER_size = 600;

btm_load = rand(size(ini_load_v));
btm_load = btm_load/sum(btm_load)*DER_size;
load_p = ini_load_p + btm_load;

v_corr_only_load_bus = loadcorrvolts(:, load_bus);
v_load_bus_conv = loadvolts(:, load_bus);
v_corr_gen_load_bus = genvolts(:, gen_load_bus);
v_load_bus = [v_corr_only_load_bus, v_corr_gen_load_bus];
v_load_bus_comp = [v_load_bus_conv, v_corr_gen_load_bus];
% Load Contingency
P_m = @(t) 0 * (t < 1) + (-cont_size/Sbase) * (t >= 1); % P_a is 0 until t = 1, then -0.1

% System Dynamics
Num = tf(1, [2*Hsys, Dsys]);
[Anum, Bnum, Cnum, Dnum] = ssdata(Num);

% Generator
Gen = tf(1, [Tg, 1]) * tf([Fh * Trh, 1], [Trh, 1]) * (1 / Rg);
% AGC = tf(5, [1, 0]);
% 
% Gen1 = parallel(Gen, AGC);
[Agen, Bgen, Cgen, Dgen] = ssdata(Gen);

% gfl
GFLsize = GFLbase + DER_size;
GFL = GFLsize/Sbase/RGFL * tf(1, [Tconv, 1]) * tf(1, [Tconv, 1])* tf(KGFL, [1, 0]) / (1+tf(KGFL, [1, 0]));
[AGFL, BGFL, CGFL, DGFL] = ssdata(GFL);

sims = 1;
for sim=1:sims
    
    % Voltage Dependency of Load    
    load_a = 0.4;%+0.05*rand();
    load_b = 0.4;%+0.05*rand();%(1-load_a).*rand(2,1);
    load_c = 1-load_a-load_b;
    loads = Sbase;
    load_riv = @(k, riv) sum(load_p.*(load_a.*(riv(k,:)./ini_load_v).^2 + load_b.*(riv(k,:)./ini_load_v) + load_c - 1))/Sbase;
    t = 0:dt:T+dt;
    
    x_gen = zeros(size(Agen, 1), N);
    x_sim = zeros(size(Anum, 1), N);
    x_gfl = zeros(size(AGFL, 1), N);
    
    y_gen = zeros(1,N);
    y_gfl = zeros(1,N);
    y_sim = zeros(1,N);

    lc = zeros(1,N);

    for i = 2:N   
        y_all = y_gen(:, i-1)+y_gfl(:, i-1)+lc(:, i-1);
        u = P_m(t(i-1))-y_all;

        [x_sim(:, i), y_sim(:, i)] = runxy(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt);
        % [x_sim(:, i), y_sim(:, i)] = runsde(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt, t(i));
        % y_sim(:, i) = y_sim(:,i) + Jump(t(i)) + Jump_rec(t(i));
        
        % volts
        lc(:,i) = load_riv(i, v_load_bus); %volts);%, volts(i,:)); %loadcorrvolts
        % lc(:,i) = load_change(t(i), y_sim(i)*60);
        [x_gen(:, i), y_gen(:, i)] = runxy(x_gen(:, i-1), y_sim(:,i), Agen, Bgen, Cgen, Dgen, dt);
        [x_gfl(:, i), y_gfl(:, i)] = runxy(x_gfl(:, i-1), y_sim(:,i), AGFL, BGFL, CGFL, DGFL, dt);
        y_gfl(:, i) = max(-0.3*GFLsize/Sbase, y_gfl(:,i));

    end

    x_gen2 = zeros(size(Agen, 1), N);
    x_sim2 = zeros(size(Anum, 1), N);
    x_gfl2 = zeros(size(AGFL, 1), N);
    
    y_gen2 = zeros(1,N);
    y_gfl2 = zeros(1,N);
    y_sim2 = zeros(1,N);

    lc2 = zeros(1,N);
    for i = 2:N   
        y_all2 = y_gen2(:, i-1)+y_gfl2(:, i-1)+lc2(:, i-1);
        u2 = P_m(t(i-1))-y_all2;

        [x_sim2(:, i), y_sim2(:, i)] = runxy(x_sim2(:, i-1), u2, Anum, Bnum, Cnum, Dnum, dt);
    
        % volts
        lc2(:,i) = load_riv(i, v_load_bus_comp);
        [x_gen2(:, i), y_gen2(:, i)] = runxy(x_gen2(:, i-1), y_sim2(:,i), Agen, Bgen, Cgen, Dgen, dt);
        [x_gfl2(:, i), y_gfl2(:, i)] = runxy(x_gfl2(:, i-1), y_sim2(:,i), AGFL, BGFL, CGFL, DGFL, dt);
        y_gfl2(:, i) = max(-0.3*GFLsize/Sbase, y_gfl2(:,i));
    end

    x_gen3 = zeros(size(Agen, 1), N);
    x_sim3 = zeros(size(Anum, 1), N);
    x_gfl3 = zeros(size(AGFL, 1), N);
    
    y_gen3 = zeros(1,N);
    y_gfl3 = zeros(1,N);
    y_sim3 = zeros(1,N);
    for i = 2:N   
        y_all3 = y_gen3(:, i-1)+y_gfl3(:, i-1);
        u3 = P_m(t(i-1))-y_all3;

        [x_sim3(:, i), y_sim3(:, i)] = runxy(x_sim3(:, i-1), u3, Anum, Bnum, Cnum, Dnum, dt);
    
        % volts
        [x_gen3(:, i), y_gen3(:, i)] = runxy(x_gen3(:, i-1), y_sim3(:,i), Agen, Bgen, Cgen, Dgen, dt);
        [x_gfl3(:, i), y_gfl3(:, i)] = runxy(x_gfl3(:, i-1), y_sim3(:,i), AGFL, BGFL, CGFL, DGFL, dt);
        y_gfl3(:, i) = max(-0.3*GFLsize/Sbase, y_gfl3(:,i));
    end
    % dvols = diff(volts);
    freq = (1+y_sim)*60;%60+y_sim*30/pi;
    freq2 = (1+y_sim2)*60;
    freq3 = (1+y_sim3)*60;
    min_freq(sim) = min(freq);

    figure(1)
    plot(t, freq, 'LineWidth', 2, 'Color', 'r', 'LineStyle', ':');
    hold on
    plot(t, freq2, 'LineWidth', 2, 'Color', 'b', 'LineStyle', '-.');
    plot(t, freq3, 'LineWidth', 2, 'Color', 'g', 'LineStyle', '--');
    if sim==sims
        plot(t, (freqsys(2:end)+1)*60, 'LineStyle', '-.', 'LineWidth', 2, 'Color', 'k')
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    fontsize(20,"points")
end

figure(2)
plot(t, -y_gen, 'LineWidth', 2);
hold on
plot(t, -y_gfl, 'LineWidth', 2);
plot(t, -lc, 'LineWidth', 2);
legend('Generator', 'GFL', 'Load')


function [x2, y2] = runxy(x, u, A, B, C, D, dt)
   drift = A*x + B*u;
   x2 = x + drift*dt;
   y2 = C*x2 + D*u;
end

function [x2, y2] = runsde(x, u, A, B, C, D, dt,t)
   drift = A*x + B*u;
   % diffusion = 0.1*u*2*(rand()-0.5)/exp(t)*(t>1);
   diffusion = 0.1*u*rand()*(t>1);
   noise = sqrt(dt)*rand(size(drift));
   x2 = x + drift*dt + diffusion*noise;
   y2 = C*x2 + D*u;
end

function sig = ini_jump_rec(dv, prev)  
    for k=1:size(prev,2)
        if (prev(1,k) == 0) && (dv(1,k)> 0.2*(10^(-3)))
            % disp(dv(1,k))
            prev(1,k)=1;
        end
    end
    sig = prev;
end