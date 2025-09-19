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

% Conventional voltage generation method based on fixed sensitivity
for i=1:size(genvolts,1)
    dloadvolts(i,:) = vsens(:, 1:29)'*dgenvolts(i,:)';
    
end
loadvolts = Volt(1, 1:29)+dloadvolts;

%% SDE Simulations
for sim=1:sims

    % BTM & Pure Load
    DER_size = 600*(1+(rand()-0.5)*0.4);
    btm_load = rand(size(ini_load_v));
    btm_load = btm_load/sum(btm_load)*DER_size; % Dirichlet distributed BTM DERs
    load_p = ini_load_p + btm_load;

    % GFL & GFM Inverters
    GFLsize = GFLbase + DER_size; % Large GFL + DER (DERs are also assumed to do GFL ctrl)
    GFMsize = 0; % Size of droop-based GFM
    VSMsize = 0; % Size of VSM GFM

    GFL = GFLsize/Sbase/RGFL * tf(1, [Tconv, 1]) * tf(1, [Tconv, 1])* tf(KGFL, [1, 0]) / (1+tf(KGFL, [1, 0]));
    GFM_droop = GFMsize/Sbase*K_GFM*tf(1, [Tconv, 1])*tf(Ome_b, [1, 0])*tf(1, [Tconv, 1])/(1+K_GFM*Droop_GFM*tf(Ome_b, [1, 0]));
    J_gfm = 2*H_gfm*VSMsize/(Ome_b^2);
    GFM_vsm = VSMsize/Sbase*K_GFM*tf(1, [Tconv, 1])*tf(Ome_b, [1, 0])*tf(1, [Tconv, 1])/(1+tf(1, [J_gfm, D_gfm])*Droop_GFM*tf(Ome_b, [1, 0]));
    
    [AGFM_droop, BGFM_droop, CGFM_droop, DGFM_droop] = ssdata(GFM_droop);
    [AGFM_vsm, BGFM_vsm, CGFM_vsm, DGFM_vsm] = ssdata(GFM_vsm);
    [AGFL, BGFL, CGFL, DGFL] = ssdata(GFL);
    
    % FFR Resources
    FFR_dev_size = 150;
    FFR_roc_size = 0;
    FFR_dev = tf(1, [Tffr_dev, 1]) * tf(1, [Tconv, 1]) * (1 / Rffr) * (FFR_dev_size)/Sbase;
    [Affrdev, Bffrdev, Cffrdev, Dffrdev] = ssdata(FFR_dev);
    
    FFR_roc = tf([1,0], 1) * tf(1, [Tffr_dev, 1]) * tf(1, [Tconv, 1]) * (1 / Rffr) * FFR_roc_size/Sbase;
    [Affrroc, Bffrroc, Cffrroc, Dffrroc] = ssdata(FFR_roc);

    % Proposed voltage generation
    vsens_rand = vsens; % Below can be used instead for voltage sensitivity with further uncertainty
    %vsens_rand = vsens-0.0003*vsens.*(rand(size(vsens)));

    dcorrload = vloadgen(H, dgenvolts, genvolts);    
    vfrat_factor = vfrat(1,1:29)/max(vfrat(1,1:29));
    v_scaling_factor = cont_size/500+0.4*(rand(29,1)-0.5)./vfrat_factor'+(cont_size/500-1)*rand(29,1);
    dcorrload2 = v_scaling_factor'.*dcorrload;

    loadcorrvolts = Volt(1, 1:29)+dcorrload2;
    v_corr_only_load_bus = loadcorrvolts(:, load_bus);
    v_corr_gen_load_bus = genvolts(:, gen_load_bus);
    v_load_bus = [v_corr_only_load_bus, v_corr_gen_load_bus];
    v_load_bus_conv = loadvolts(:, load_bus);
    v_load_bus_comp = [v_load_bus_conv, v_corr_gen_load_bus];

    % Voltage Dependency of Load (ZIP Load Model)
    load_a = 0.4-0.2*(rand()-0.5); % Coefficient Z
    load_b = 0.4-0.2*(rand()-0.5); % Coefficient I
    load_c = 1-load_a-load_b;      % Coefficient P
    load_riv = @(k, riv) sum(load_p.*(load_a.*(riv(k,:)./ini_load_v).^2 + load_b.*(riv(k,:)./ini_load_v) + load_c - 1))/(sum(load_p)+cont_size);
    
    t = 0:dt:T+dt;
    
    x_gen = zeros(size(Agen, 1), N);
    x_sim = zeros(size(Anum, 1), N);
    x_gfl = zeros(size(AGFL, 1), N);
    x_gfm = zeros(size(AGFM_droop, 1), N);
    x_vsm = zeros(size(AGFM_vsm, 1), N);
    x_ffrdev = zeros(size(Affrdev, 1), N);
    x_ffrroc = zeros(size(Affrroc, 1), N);
    trip_size = zeros(1,N);

    y_gen = zeros(1,N);
    y_gfl = zeros(1,N);
    y_gfm = zeros(1,N);
    y_vsm = zeros(1,N);
    y_sim = zeros(1,N);
    y_ffrdev = zeros(1,N);
    y_ffrroc = zeros(1,N);

    lc = zeros(1,N);
    lc2 = zeros(1,N);
    for i = 2:N   
        % 
        y_all = y_gen(:, i-1)+y_gfl(:, i-1)+lc(:, i-1)+y_gfm(:, i-1)+y_vsm(:, i-1)+y_ffrdev(:, i-1)+y_ffrroc(:, i-1);
        u = P_m(t(i-1))-y_all;
        
        % To include unknown uncertainty, use runsde; Else, use runxy
        [x_sim(:, i), y_sim(:, i)] = runxy(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt);
        % [x_sim(:, i), y_sim(:, i)] = runsde(x_sim(:, i-1), u, Anum, Bnum, Cnum, Dnum, dt, t(i));

        % Can be added to include jump uncertainty
        % y_sim(:, i) = y_sim(:,i) + Jump(t(i)) + Jump_rec(t(i));
        
        % Voltages & DER Trips
        v_load_bus(i,:) = v_load_bus(i,:)-0.1*v_load_bus(i,:).*rand(1,19).*(trip_size(i));
        lc(:,i) = load_riv(i, v_load_bus);
        trip_signal(i,:) = DER_vrt(size(v_load_bus,2), t(i), dt, v_load_bus(i,:)); % All DERs are assumed to be legacy DER in here
        trip_size(i) = sum(btm_load.*(trip_signal(i,:)>0));

        [x_gen(:, i), y_gen(:, i)] = runxy(x_gen(:, i-1), y_sim(:,i), Agen, Bgen, Cgen, Dgen, dt);
        [x_gfl(:, i), y_gfl(:, i)] = runxy(x_gfl(:, i-1), y_sim(:,i), AGFL, BGFL, CGFL, DGFL, dt);
        [x_gfm(:, i), y_gfm(:, i)] = runxy(x_gfm(:, i-1), y_sim(:,i), AGFM_droop, BGFM_droop, CGFM_droop, DGFM_droop, dt);
        [x_vsm(:, i), y_vsm(:, i)] = runxy(x_vsm(:, i-1), y_sim(:,i), AGFM_vsm, BGFM_vsm, CGFM_vsm, DGFM_vsm, dt);
        [x_ffrdev(:,i), y_ffrdev(:,i)] = runxy(x_ffrdev(:, i-1), y_sim(:,i), Affrdev, Bffrdev, Cffrdev, Dffrdev, dt);
        [x_ffrroc(:,i), y_ffrroc(:,i)] = runxy(x_ffrroc(:, i-1), y_sim(:,i), Affrroc, Bffrroc, Cffrroc, Dffrroc, dt);

        y_gfl(:, i) = max(-0.2*GFLsize/Sbase, y_gfl(:,i));
        y_gfl(:, i) = y_gfl(:,i)*(GFLsize-trip_size(i))/GFLsize;
        y_gfm(:, i) = max(-0.2*GFMsize/Sbase, y_gfm(:,i));
        y_vsm(:, i) = max(-0.2*VSMsize/Sbase, y_vsm(:,i));

    end
    freq = (1+y_sim)*60;
    min_freq(sim) = min(freq);

    figure(1)
    plot(t, freq, 'LineStyle', ':', 'LineWidth', 2, 'Color', 'r');
    hold on
    if sim==sims
        plot(t, (freqsys(2:end)+1)*60, 'Color', 'k', 'LineWidth', 2)
        hold on
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    fontsize(24,"points")
end

% S_LG and S_GG figures
figure(2)
heatmap(vsens(:,1:29)')
colormap('jet')
xlabel('Generator bus')
ylabel('Load bus')
title('Generator to Load Sensitivity')

figure(3)
heatmap(vsens(:,30:39)')
colormap('jet')
xlabel('Generator bus')
ylabel('Generator bus')
title('Generator to Generator Sensitivity')

figure(4)
subplot(2,1,1)
plot(time, (Freq+1)*60, 'LineWidth', 2)
hold on 
plot(time, (freqsys+1)*60, 'LineWidth', 3, 'Color', 'k')
ylabel('Frequency (Hz)');
title('PSSE Results')
subplot(2,1,2)
plot(time, Volt, 'LineWidth', 2)
xlabel('Time (sec)')
ylabel('Voltage (p.u.)');
fontsize(24, 'points')

figure(5)
plot(t, -y_gen, 'LineWidth', 3);
hold on
plot(t, -y_gfl, 'LineWidth', 3, 'LineStyle', '--');
plot(t, -lc, 'LineWidth', 3, 'LineStyle', '-.');
legend('Generator', 'GFL', 'Load')
xlabel('Time (s)');
ylabel('Normalized P');
fontsize(24,"points")

figure(6)
for i=1:size(load_bus,2)
    a = load_bus(i);
    plot(time, Volt(:,a), 'LineWidth', 2, 'Color','k')
    hold on
    plot(time, v_load_bus(:,i), 'r--', 'LineWidth', 2)
    plot(time, v_load_bus_conv(:,i), 'b:', 'LineWidth', 2)
end


[mu,s,muci,sci] = normfit(min_freq);
% gamma = fitdist(min_freq', 'Gamma');

x = min(min_freq):0.001:max(min_freq);
y = normpdf(x,mu,s);
% y_gamma = gampdf(x,gamma.a,gamma.b);

pd = fitdist(min_freq', 'Normal');
ci = paramci(pd,'Alpha',.01);   

figure(7)
h = histogram(min_freq, 30, 'Orientation', 'horizontal','Normalization','pdf');
hold on
plot(y,x, 'LineWidth', 2)
yline(min((freqsys(2:end)+1)*60), 'r')
% hold on
% plot(x,y_gamma, 'LineWidth', 2, 'LineStyle', '--')
ylabel('Frequency Nadir');
fontsize(20,"points")

figure(8)
h = histogram(min_freq, 30, 'Normalization','pdf');
hold on
plot(x,y, 'LineWidth', 2)
%yline(min((freqsys(2:end)+1)*60), 'r')
% hold on
% plot(x,y_gamma, 'LineWidth', 2, 'LineStyle', '--')
ylabel('Frequency Nadir');
fontsize(20,"points")

disp([mu-3*s, mu-2*s])

function dcorrload=vloadgen(H, dgenvolts, genvolts)
    H_G = H(30:end, :);
    [u, s, v] = svd(H);
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
        dcorrload(i,:) = S_LGt*dV';
    end
end

function [x2, y2] = runxy(x, u, A, B, C, D, dt)
   drift = A*x + B*u;
   x2 = x + drift*dt; 
   y2 = C*x2 + D*u;
end

function [x2, y2] = runsde(x, u, A, B, C, D, dt,t)
   drift = A*x + B*u;
   % diffusion = 0.1*u*2*(rand()-0.5)/exp(t)*(t>1);
   diffusion_f = (0.003*(rand()-0.5)+0.3*x*(rand()-0.5))*(t>1);
   noise_f = sqrt(dt)*rand(size(drift));
   diffusion_p = (0.0015*(rand()-0.5)+0.15*u*(rand()-0.5))*(t>1);
   noise_p = sqrt(dt)*rand(size(drift));
   x2 = x + drift*dt+ diffusion_p*noise_p;
   y2 = C*x2 + D*u + diffusion_f*noise_f;
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

function vrt_out = DER_vrt(N, t, dt, loadvolt)
    persistent v_buffer
    persistent der_trip_rec

    buffer_size_1 = round(0.5/dt); % 0.5 sec under 0.88 p.u.

    if isempty(v_buffer)
        v_buffer = zeros(N,1);
    end

    if isempty(der_trip_rec)
        der_trip_rec = zeros(N,1);
    end
    if t<0.1
        v_buffer = zeros(N,1);
        der_trip_rec = zeros(N,1);
    end

    for b=1:size(loadvolt, 2)
        if der_trip_rec(b)==1
            der_trip_rec(b)=1;
        else
            uv1 = (loadvolt(b)<0.88);

            if uv1
                v_buffer(b,1) = v_buffer(b,1)+1;
            else
                v_buffer(b,:) = 0;
            end

            if (v_buffer(b,1)>buffer_size_1)
                der_trip_rec(b) = 1;
            end
        end
    end
    vrt_out = der_trip_rec;
end
