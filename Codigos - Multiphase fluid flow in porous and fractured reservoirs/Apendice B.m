%B.1 MATLAB CODE OF WELGE GRAPHIC METHOD IN A RADIAL SYSTEM

%% Problem description
% Analytical solution of Buckley–Leverett equation
% in a radial system using Welge method
% water phase fluid displacing oil phase fluid
% Brooks–Corey type relative permeability curves are used

close all; clear all; clc;
%% 1 Parameters initialization
% 1.1 rock properties
R = 10.00; % length of domain [m]
r = 0.10; % radii of injection well [m]
h = 1.0; % height of domain [m]
phi = 0.25; % porosity
k = 9.869e-14; % absolute permeability [m^2]

% 1.2 fluid properties
muo = 5.0e-3; % oil phase viscosity [Pa*s]
muw = 1.0e-3; % water phase viscosity [Pa*s]
rho_o = 0.8e3; % oil density [kg/m^3]
rho_w = 1.0e3; % water density [kg/m^3]
dlt_rho = rho_w - rho_o; % density difference [kg/m^3]

% 1.3 relative permeability parameters
Sor = 0.15; % residual oil saturation
Swc = 0.15; % connate water saturation
no = 2.00; % exponent of oil phase
nw = 2.00; % exponent of water phase
kro_max = 0.75; % maximum permeability of oil
krw_max = 0.75; % maximum permeability of water

% 1.4 other initial parameters
dlt_Sw = 1e-3; % constant saturation step
Sw = [(Swc + eps):dlt_Sw:(1 - Sor)]'; % water saturation vector
qt = 1.0e-4; % constant water injection rate [m^3/s]
time0 = 86400 * 1.0; % initial calculation time [s]
timef = 86400 * 2.2796; % final calculation time [s]
nts = 4; % time steps for calculating
nsw = size(Sw, 1); % number of water saturation vector
format = '%4.2e'; % precision format for plotting legend

%% 2 Relative permeabilities, fractional flow function and its derivative
% 2.1 initialization of function handles
kr_w = @(S) krw_max .* ((S - Swc) / (1 - Swc - Sor)) .^ nw; % relative permeability of water
kr_o = @(S) kro_max .* ((1 - S - Sor) / (1 - Swc - Sor)) .^ no; % relative permeability of oil
mob = @(S) (kr_o(S) ./ muo) .* (muw ./ kr_w(S)); % mobility function
f_w = @(S) 1 ./ (1 + mob(S)); % fractional flow function fw
df_w = @(S) ((nw .* mob(S) ./ (S - Swc)) + (no .* mob(S)) ./ (1 - S - Sor))... % derivative of fw
./ (1 + mob(S)) .^ 2;

% 2.2 relative permeability and fractional flow vectors
krw = kr_w(Sw); % relative permeability vector of water
kro = kr_o(Sw); % relative permeability vector of oil
fw = f_w(Sw); % fractional flow vector
dfw = df_w(Sw); % fractional flow derivate vector

%% 3 Calculate advance frontal water saturation
% 3.1 Evaluate advance frontal water saturation Swf
[index, Swf, dfwf, bf] = calSatFront(Sw, fw, dfw);
% 3.2 Calculate the time when Swf reaches at production well
krw_swf = kr_w(Swf);
kro_swf = kr_o(Swf);
fw_swf = f_w(Swf);
dfw_swf = df_w(Swf);
t_pw = pi * h * phi * (R ^ 2 - r ^ 2) / (qt * dfw_swf);

%% 4 Calculate water saturation profile
t = linspace(time0, timef, nts);
dfwt = dfw(end:-1:index); % inverted sequence of dfw
Swt = Sw(end:-1:index); % inverted sequence of Sw

for ti = 1:nts

    for i = 1:(nsw - index + 1)
        Rsw(i, ti) = sqrt(r ^ 2 + qt * t(ti) / (pi * h * phi) * dfwt(i));
    end

end

%% 5 Plot results
% 5.2 Plot the relative permeability curves
h_fig1 = figure(1);
set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves:kr');
plot(Sw, krw, '-b', Sw, kro, '-r', 'LineWidth', 2.0);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Brooks–Corey relative permeability curves')
xlabel('{\it S}_w');
ylabel('{\it k}_{r {\it \beta}}');
set(gca, 'YTick', 0:0.2:1);
h_legend1 = legend('water phase', 'oil phase');
set(h_legend1, 'Box', 'on', 'Location', 'best');

% 5.2 Plot the fractional flow function curve
h_fig2 = figure(2);
set(h_fig2, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional Flow Curve fw');
plot(Sw, fw, '-b', 'LineWidth', 2.0);
hold on;
plot(Sw, (dfwf .* Sw + bf), '-r', 'LineWidth', 2.0);
plot(Swc, 0, 'ro', 'Markersize', 8);
plot(Swf, fw(index), 'ro', 'Markersize', 8);
plot([Swf Swf], [0 fw(index)], '--r', 'LineWidth', 1.5);
plot([0 Swf], [fw(index) fw(index)], '--r', 'LineWidth', 1.5);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('fractional flow function curve')
xlabel('{\it S}_w');
ylabel('{\it f}_w');
set(gca, 'YTick', 0:0.2:1);

% 5.3 Plot the derivative curve of fractional flow function
h_fig3 = figure(3);
set(h_fig3, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Derivatives of Fractional Flow Curve dfw / dSw');
plot(Sw, dfw, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('derivatives of fractional flow function')
xlabel('{\it S}_w');
ylabel('d {\it f}_w / d {\it S}_w');
axis square;

% 5.4 Plot the water saturation profiles
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'water Saturation Profiles:Sw(t)');
SwTime = [];

for ti = 1:nts
    plot(Rsw(:, ti), Swt, '-r', 'LineWidth', 2.0);
    hold on;
    plot([0 Rsw(1, ti)], [1 - Sor 1 - Sor], '-r', 'LineWidth', 2.0);
    SwTime = [SwTime; ['Time = ' num2str(t(1, ti) / 86400, format) ' day']];
    plot([Rsw(end, ti) Rsw(end, ti)], [Swf Swc], '-r', 'LineWidth', 2.0);
    plot([Rsw(end, ti) R], [Swc Swc], '-r', 'LineWidth', 2.0);
end

set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 R]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('water saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_w');
h_legend2 = legend(SwTime);
set(h_legend2, 'Box', 'on', 'Location', 'best');
axis square;
%% end

% B.2 MATLAB CODE OF INTEGRAL METHOD BASED ON
% MASS BALANCE PRINCIPLE IN A RADIAL SYSTEM
%% Problem description
% Analytical solution of Buckley–Leverett equation
% in a radial system using integral method
% wetting phase fluid displacing non-wetting phase fluid
% Brooks–Corey type relative permeability curves are used
close all; clear all; clc;

%% 1 Parameters initialization
% 1.1 rock properties
R = 100.0; % length of domain [m]
r = 0.02; % radii of injection well [m]
h = 1.0; % height of domain [m]
phi = 0.25; % porosity
k = 9.869e-13; % absolute permeability [m^2]

% 1.2 fluid properties
muo = 5.0e-3; % oil phase viscosity [Pa*s]
muw = 1.0e-3; % water phase viscosity [Pa*s]
rho_o = 0.8e3; % oil density [kg/m^3]
rho_w = 1.0e3; % water density [kg/m^3]
dlt_rho = rho_w - rho_o; % density difference [kg/m^3]

% 1.3 relative permeability parameters
Sor = 0.15; % residual oil saturation
Swc = 0.15; % connate water saturation
no = 1.00; % exponent of oil phase
nw = 2.00; % exponent of water phase
kro_max = 0.85; % maximum permeability of oil phase
krw_max = 0.85; % maximum permeability of water phase

% 1.4 other initial parameters
dlt_Sw = 1e-3; % constant saturation step
Sw = [(Swc + eps):dlt_Sw:(1 - Sor)]'; % water saturation vector
qt = 2.0e-3; % constant water injection rate [m^3/s]
time0 = 86400 * 1.0; % initial calculation time
timef = 86400 * 9.0; % final calculation time
nts = 4; % time steps for calculating
nsw = size(Sw, 1); % number of wetting-phase saturation
format = '%3.2f'; % precision format for plotting legend

%% 2 Relative permeabilities, fractional flow function and its derivative
% 2.1 initialization of function handles
kr_w = @(S) krw_max .* ((S - Swc) / (1 - Swc - Sor)) .^ nw; % relative permeability of water
kr_o = @(S) kro_max .* ((1 - S - Sor) / (1 - Swc - Sor)) .^ no; % relative permeability of oil
mob = @(S) (kr_o(S) ./ muo) .* (muw ./ kr_w(S)); % mobility function
f_w = @(S) 1 ./ (1 + mob(S)); % fractional flow function fw
df_w = @(S) ((nw .* mob(S) ./ (S - Swc)) + (no .* mob(S)) ./ (1 - S - Sor))... % derivative of fw
./ (1 + mob(S)) .^ 2;

% 2.2 relative permeability and fractional flow vectors
krw = kr_w(Sw); % relative permeability vector of water
kro = kr_o(Sw); % relative permeability vector of oil
fw = f_w(Sw); % fractional flow vector
dfw = df_w(Sw); % fractional flow derivate vector

%% 3 Calculate advance frontal saturation and travelling distance
% 3.1 calculate travelling distance of every saturations
t = linspace(time0, timef, nts);
dfwt = dfw(end:-1:1); % inverted sequence of dfw
fwt = fw(end:-1:1); % inverted sequence of fw
Swt = Sw(end:-1:1); % inverted sequence of Sw

for ti = 1:nts
    Wi(ti) = qt * t(ti); % injected wetting phase fluid volume

    for i = 1:nsw
        Rsw(i, ti) = sqrt(r ^ 2 + qt * t(ti) / (pi * h * phi) * dfwt(i)); % the travelling distance
    end

    Rsw0 = [r; Rsw(1:end - 1, ti)]; % the 2nd travelling distance vector
    dlt_Rsw = Rsw(:, ti) .^ 2 - Rsw0 .^ 2; % dlt_Rsw=Rsw_j^2 - Rsw_(j-1)^2, r_0=r
    dlt_Swt = Swt - Swc; % dlt_Swt = Swt,j - Swc

    for i = 1:nsw
        % calculate the injected fluid volume from r to Rswf
        Vi = pi * phi * h * sum(dlt_Rsw(1:i) .* dlt_Swt(1:i));

        if Vi >= Wi(ti)
            index(ti) = i; % index of the Swf in Swt
            Swf (ti) = Swt(i); % advance front saturation: Swf
            Rswf (ti) = Rsw(i, ti); % the position of Swf
            break;
        end

    end

end

Rsw = Rsw(1:index(1), :);
Swt = Swt(1:index(1));

% 3.2 Calculate the time when Swf reaches at production well
krw_swf = kr_w(Swf(1));
kro_swf = kr_o(Swf(1));
fw_swf = f_w(Swf(1));
dfw_swf = df_w(Swf(1));
t_pw = pi * h * phi * (R ^ 2 - r ^ 2) / (qt * dfw_swf);

%% 4 Plot results
% 4.1 Plot the relative permeability curves
h_fig1 = figure(1);
set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves:kr');
plot(Sw, krw, '-b', Sw, kro, '-r', 'LineWidth', 2.0);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Brooks–Corey relative permeability curves')
xlabel('{\it S}_w');
ylabel('{\it k}_{r {\it \beta}}');
set(gca, 'YTick', 0:0.2:1);
h_legend1 = legend('water phase', 'oil phase');
set(h_legend1, 'Box', 'on', 'Location', 'best');

% 4.2 Plot the fractional flow function curve
h_fig2 = figure(2);
set(h_fig2, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional Flow Curve fw');
plot(Sw, fw, '-b', 'LineWidth', 2.0);
hold on;
b = fw_swf - dfw_swf * Swf(1);
plot(Sw, (dfw_swf .* Sw + b), '-r', 'LineWidth', 2.0);
plot(Swc, 0, 'ro', 'Markersize', 8);
plot(Swf(1), fw_swf, 'ro', 'Markersize', 8);
plot([Swf(1) Swf(1)], [0 fwt(index(1))], '--r', 'LineWidth', 1.5);
plot([0 Swf(1)], [fwt(index(1)) fwt(index(1))], '--r', 'LineWidth', 1.5);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('fractional flow function curve')
xlabel('{\it S}_w');
ylabel('{\it f}_w');
set(gca, 'YTick', 0:0.2:1);

% 4.3 Plot the derivative curve of fractional flow function
h_fig3 = figure(3);
set(h_fig3, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Derivatives of Fractional Flow Curve dfw / dSw');
plot(Sw, dfw, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('derivatives of fractional flow function')
xlabel('{\it S}_w');
ylabel('d {\it f}_w / d {\it S}_w');
axis square;

% 4.4 Plot the water saturation profiles
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'water Saturation Profiles: Sw(t)');
SwTime = [];

for ti = 1:nts
    plot(Rsw(:, ti), Swt, '-r', 'LineWidth', 2.0);
    hold on;
    plot([r Rsw(1, ti)], [1 - Sor 1 - Sor], '-r', 'LineWidth', 2.0);
    SwTime = [SwTime; ['Time = ' num2str(t(1, ti) / 86400, format) ' day']];
    plot([Rsw(end, ti) Rsw(end, ti)], [Swf(ti) Swc], '-r', 'LineWidth', 2.0);
    plot([Rsw(end, ti) R], [Swc Swc], '-r', 'LineWidth', 2.0);
end

set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
set(gca, 'YLim', [0 1]);
title('water saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_w');
h_legend4 = legend(SwTime);
set(h_legend4, 'Box', 'on', 'Location', 'best');
axis square;
%% end

% B.3 MATLAB CODE FOR BUCKLEY–LEVERETT SOLUTION
% IN A LINEAR COMPOSITE SYSTEM
%% Problem description
% Analytical solution of Buckley–Leverett equation using mass balance method
% in one-dimensional composite porous medium system
% Brooks–Corey type relative permeability curves are used
close all; clear all; clc;
%% 1 Parameters initialization for two domains
% 1.1 rock properties in domain 1
L1 = 6.00; % length of domain 1 [m]
phi1 = 0.30; % porosity of domain 1
k1 = 1.0e-14; % permeability of domain 1 [m^2]
Sor1 = 0.20; % residual oil saturation in domain 1
Swc1 = 0.20; % connate water saturation in domain 1
no1 = 2.50; % exponent of oil phase in domain 1
nw1 = 1.50; % exponent of water phase in domain 1
kro1_max = 0.80; % maximum permeability of oil domain 1
krw1_max = 0.80; % maximum permeability of water domain 1
dlt_Sw1 = 1e-3; % constant saturation step in domain 1
Sw1 = [(Swc1 + eps):dlt_Sw1:(1 - Sor1 - eps)]'; % water saturation vector for domain 1
% 1.2 rock properties in domain 2
L2 = 6.00; % length of domain 2 [m]
phi2 = 0.30; % porosity of domain 2
k2 = 1.0e-14; % permeability of domain 2 [m^2]
Sor2 = 0.20; % residual oil saturation in domain 2
Swc2 = 0.20; % connate water saturation in domain 2
no2 = 1.50; % exponent of oil phase in domain 2
nw2 = 2.50; % exponent of water phase in domain 2
kro2_max = 0.75; % maximum permeability of oil domain 2
krw2_max = 0.75; % maximum permeability of water domain 2
dlt_Sw2 = 1e-3; % constant saturation step in domain 2
Sw2 = [(Swc2 + eps):dlt_Sw2:(1 - Sor2 - eps)]'; % water saturation vector for domain 2

    % 1.3 fluid properties for the composite system
    muo = 5.0e-3; % oil phase viscosity [Pa*s]
    muw = 1.0e-3; % water phase viscosity [Pa*s]
    rho_o = 0.8e3; % oil density [kg/m^3]
    rho_w = 1.0e3; % water density [kg/m^3]
    dlt_rho = rho_w - rho_o; % density difference [kg/m^3]
    % 1.4 other parameters for calculation
    A = 1.0; % area of cross-section [m^2]
    qt = 1.0e-5; % constant water injection rate [m^3/s]
    g = 9.8067; % gravity acceleration constant [m/s^2]
    alpha = pi * 0.0; % angle between x and horizontal [rad]
    time0 = 86400 * 0.1; % initial calculation time [s]
    timef = 86400 * 1.0; % final calculation time [s]
    nts1 = 5; % time steps for domain 1
    nts2 = 5; % time steps for domain 2 [-]
    nsw1 = size(Sw1, 1); % number of water saturation in domain 1
    nsw2 = size(Sw2, 1); % number of water saturation in domain 2
    format = '%4.2e'; % precision format for plotting
    %% 2 Calculate fractional flow function and plots
    % 2.1 initialization of function handles in domain 1
    kr_w1 = @(S) krw1_max .* ((S - Swc1) / (1 - Swc1 - Sor1)) .^ nw1; % relative permeability krw1
    kr_o1 = @(S) kro1_max .* ((1 - S - Sor1) / (1 - Swc1 - Sor1)) .^ no1; % relative permeability kro1
    mob1 = @(S) (kr_o1(S) ./ muo) .* (muw ./ kr_w1(S)); % mobility function
    f_w1 = @(S) 1 ./ (1 + mob1(S)) - A * k1 .* kr_o1(S) ./ (qt * muo) * ... % fractional flow fw1
    dlt_rho * g * sin(alpha) ./ (1 + mob1(S));
    df_w1 = @(S) ((nw1 .* mob1(S) ./ (S - Swc1)) + (no1 .* mob1(S)) ./ (1 - S - Sor1))... % derivate, fw1
    ./ (1 + mob1(S)) .^ 2 + A * k1 .* kr_o1(S) * no1 * dlt_rho * g * sin(alpha) ./ ...
    ((1 - S - Sor1) * qt * muo .* (1 + mob1(S))) + A * k1 .* kr_o1(S) * dlt_rho * g * sin(alpha) ./ ...
    (qt * muo .* (1 + mob1(S)) .^ 2) .* (((nw1 .* mob1(S) ./ (S - Swc1)) + (no1 .* mob1(S))...
    ./ (1 - S - Sor1)) ./ (1 + mob1(S)) .^ 2);
    % 2.1 initialization of function handles in domain 2
    kr_w2 = @(S) krw2_max .* ((S - Swc2) / (1 - Swc2 - Sor2)) .^ nw2; % relative permeability krw2
    kr_o2 = @(S) kro2_max .* ((1 - S - Sor2) / (1 - Swc2 - Sor2)) .^ no2; % relative permeability kro2
    mob2 = @(S) (kr_o2(S) ./ muo) .* (muw ./ kr_w2(S)); % mobility function
    f_w2 = @(S) 1 ./ (1 + mob2(S)) - A * k2 .* kr_o2(S) ./ (qt * muo) * ... % fractional flow fw2
    dlt_rho * g * sin(alpha) ./ (1 + mob2(S));
    df_w2 = @(S) ((nw2 .* mob2(S) ./ (S - Swc2)) + (no2 .* mob2(S)) ./ (1 - S - Sor2))... % derivate fw2
    ./ (1 + mob2(S)) .^ 2 + A * k2 .* kr_o2(S) * no2 * dlt_rho * g * sin(alpha) ./ ...
    ((1 - S - Sor2) * qt * muo .* (1 + mob2(S))) + A * k2 .* kr_o2(S) * dlt_rho * g * sin(alpha) ./ ...
    (qt * muo .* (1 + mob2(S)) .^ 2) .* (((nw2 .* mob2(S) ./ (S - Swc2)) + (no2 .* mob2(S))...
    ./ (1 - S - Sor2)) ./ (1 + mob2(S)) .^ 2);
    % 2.3 relative permeability and fractional flow vectors in domain 1
    krw1 = kr_w1(Sw1); % relative permeability vector of water
    kro1 = kr_o1(Sw1); % relative permeability vector of oil
    fw1 = f_w1(Sw1); % fractional flow vector
    dfw1 = df_w1(Sw1); % fractional flow derivate vector
    % 2.4 relative permeability and fractional flow vectors in domain 2
    krw2 = kr_w2(Sw2); % relative permeability vector of water
    kro2 = kr_o2(Sw2); % relative permeability vector of oil
    fw2 = f_w2(Sw2); % fractional flow vector
    dfw2 = df_w2(Sw2); % fractional flow derivate vector
    % 2.5 plot relative permeability curves
    h_fig1 = figure(1);
    set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves: kr');
    plot(Sw1, krw1, '-b', Sw1, kro1, '--b', 'LineWidth', 2.0);
    hold on;
    plot(Sw2, krw2, '-r', Sw2, kro2, '--r', 'LineWidth', 2.0);
    axis([0.0 1.0 0.0 1.0]);
    axis square;
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('Brooks–Corey relative permeability curves')
    xlabel('{\it S}_w');
    ylabel('{\it k}_{r{\it \beta}}');
    set(gca, 'YTick', 0:0.2:1);
    h_legend1 = legend('domain 1: water', 'domain 1: oil', 'domain 2: water', 'domain 2: oil');
    set(h_legend1, 'Box', 'on', 'Location', 'best');
    %% 3 Calculate advance frontal water saturation in domain 1 and plot
    % 3.1 calculate sock front of water saturation
    [index, Swf1, dfwf1, bf1] = calSatFront(Sw1, fw1, dfw1);
    % 3.2 Plot the fractional flow function curves
    h_fig2 = figure(2);
    set(h_fig2, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional Flow Curves: fw');
    plot(Sw1, fw1, '-b', 'LineWidth', 2.0);
    hold on;
    plot(Sw2, fw2, '--b', 'LineWidth', 2.0);
    plot(Sw1, (dfwf1 .* Sw1 + bf1), '-r', 'LineWidth', 2.0);
    plot(Swc1, 0, 'ro');
    plot(Swf1, fw1(index), 'ro');
    plot([Swf1 Swf1], [0 fw1(index)], '--r', 'LineWidth', 2.0);
    plot([0 Swf1], [fw1(index) fw1(index)], '--r', 'LineWidth', 2.0);
    axis([0.0 1.0 0.0 1.0]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('fractional flow function curves');
    xlabel('{\it S}_w');
    ylabel('{\it f}_w');
    h_legend2 = legend('domain 1', 'domain 2');
    set(gca, 'YTick', 0:0.2:1);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    set(h_legend2, 'Box', 'on', 'Location', 'best');
    axis square;
    % 3.3 Plot the derivative curve of fractional flow function
    h_fig3 = figure(3);
    set(h_fig3, 'color', 'w', 'NumberTitle', 'off','Name', 'Derivatives of Fractional Flow Curves: dfw/dSw');
    plot(Sw1, dfw1, '-b', 'LineWidth', 2.0);
    hold on;
    plot(Sw2, dfw2, '--b', 'LineWidth', 2.0);
    set(gca, 'XLim', [0 1]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('derivatives of fractional flow function');
    xlabel('{\it S}_w');
    ylabel('d {\it f}_w / d {\it S}_w');
    h_legend3 = legend('domain 1', 'domain 2');
    set(h_legend3, 'Box', 'on', 'Location', 'best');
    axis square;
    %% 4 Water saturation profiles in domain 1 before Swf1 reaches interface
    % 4.1 the time when Swf reaches at the interface between domain1 and domain 2
    krw1_swf = kr_w1(Swf1);
    kro1_swf = kr_o1(Swf1);
    fw1_swf = f_w1(Swf1);
    dfw1_swf = df_w1(Swf1);
    t_inf = A * phi1 * L1 / (qt * dfw1_swf);
    % 4.2 calculate the travelling distance of specific water saturations Sw1
    t1 = linspace(time0, t_inf, nts1); % time vector for domain 1
    dfwt1 = dfw1(end:-1:index); % inverted sequence of dfw1
    Swt1 = Sw1(end:-1:index); % inverted sequence of Sw1

    for ti = 1:nts1 % loop time vector

        for i = 1:(nsw1 - index + 1)
            Xsw1(i, ti) = qt * t1(ti) / (A * phi1) * dfwt1(i);
        end

    end

    % 4.3 figure 4: Plot the water saturation profiles in domain 1
    h_fig4 = figure(4);
    set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Water Saturation Profiles');
    SwTime = [];

    for ti = 1:nts1
        plot(Xsw1(:, ti), Swt1, '-r', 'LineWidth', 2.0);
        hold on;
        plot([0 Xsw1(1, ti)], [1 - Sor1 1 - Sor1], '-r', 'LineWidth', 2.0);
        SwTime = [SwTime; ['Time = ' num2str(t1(ti) / 86400, format) ' day']];
        plot([Xsw1(end, ti) Xsw1(end, ti)], [Swf1 Swc1], '-r', 'LineWidth', 2.0);
        plot([Xsw1(end, ti) L1], [Swc1 Swc1], '-r', 'LineWidth', 2.0);
    end

    set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L1]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('water saturation profiles in domain 1');
    xlabel('{\it x} (m)'); ylabel('{\it S}_w');
    h_legend4 = legend(SwTime);
    set(h_legend4, 'Box', 'on', 'Location', 'best');
    axis square;

    %% 5 Calculate water saturation profile in domain 1 when time > t_inf
    % 5.1 calculate the corresponding parameters at interface
    if timef <= t_inf
        error('Please reset the timef due to timef <= t_inf!');
    end

    t2 = linspace(t_inf, timef, nts2)';
    dfw1_inf = A * phi1 * L1 / qt ./ t2; % dfw1 at the interface (upstream)
    Sw1_inf = calObjFun(dfw1_inf, Sw1, dfw1, index); % Sw,1 at the interface(upstream)
    krw1_inf = kr_w1(Sw1_inf); % relative permeability of water
    kro1_inf = kr_o1(Sw1_inf); % relative permeability of oil
    fw1_inf = f_w1(Sw1_inf); % fractional flow function
    Sw1_avg = Sw1_inf + (1 - fw1_inf) ./ dfw1_inf; % average saturation of domain 1
    % 5.2 % calculate the travelling distance of water saturations
    dfwt1 = dfw1(end:-1:index); % inverted sequence of dfw1
    Swt1 = Sw1(end:-1:index); % inverted sequence of Sw1

    for ti = 1:nts2 % loop times vector

        for i = 1:(nsw1 - index + 1)
            Xsw1(i, ti + nts1) = qt * t2(ti) / (A * phi1) * dfwt1(i);
        end

    end

    %% 6 Calculate water saturation profile in domain 2 when time > t_inf
    % 6.1 calculate water saturation at interface and injected fluid
    fw2_inf = fw1_inf; % fw,2 at the interface (downstream)
    Sw2_inf = calObjFun(fw2_inf, Sw2, fw2); % Sw,2 at the interface (downstream)
    Qt = qt .* t2; % the injected water into domains
    W1 = A * phi1 * L1 .* (Sw1_avg - Swc1); % the injected water into domain 1
    W2 = Qt - W1; % the injected water into domain 2
    if (abs(W2(1) / Qt(1)) <= 1e-3); W2(1) = 0; end; % enforce W2 == 0.0
    % 6.1 calculate saturation profiles using mass balance principle
    Sw2_k = cell(nts2, 1); % water saturation Sw,2: initialization
    Xsw2_k = cell(nts2, 1); % initial travelling distance of Sw,2
    Sw2_k{1, 1} = [Sw2_inf(1); Swc2]; % the saturation Sw,2 set at t_inf
    Xsw2_k{1, 1} = [L1; L1]; % the distance Xsw,2 set at t_inf

    for ti = 2:nts2
        % 6.1.1 select a set of saturation Sw,2
        % choose a set of saturation Sw,2: [Sw2,inf(1) Sw2,inf(ti)]
        Sw2_k1 = [(Sw2_inf(ti) - dlt_Sw2):(-dlt_Sw2):Sw2_inf(1)]';
        % choose a set of saturation Sw,2: [Swc2 Sw2,inf(1)])
        Sw2_k2 = [(Sw2_inf(1) - dlt_Sw2):(-dlt_Sw2):Swc2]';
        % the whole set of water saturation Sw,2
        Sw2_k{ti, 1} = [Sw2_inf(ti); Sw2_k1; Sw2_k2];
        % 6.1.2 calculate saturation profiles for Sw2,inf(1) <= Sw,2 <= Sw2,inf(ti)
        krw2_k1 = kr_w2(Sw2_k1); % relative permeability of water
        kro2_k1 = kr_o2(Sw2_k1); % relative permeability of oil
        fw2_k1 = f_w2(Sw2_k1); % fractional flow function
        dfw2_k1 = df_w2(Sw2_k1); % derivate of fractional flow function
        fw1_k1 = fw2_k1;
        Sw1_k1 = calObjFun(fw1_k1, Sw1, fw1); % Sw1 at interface in domain 1
        krw1_k1 = kr_w1(Sw1_k1); % relative permeability of water
        kro1_k1 = kr_o1(Sw1_k1); % relative permeability of oil
        fw1_k1 = f_w1(Sw1_k1); % fractional flow function
        dfw1_k1 = df_w1(Sw1_k1); % derivate of fractional flow function
        ts = A * phi1 * L1 / qt ./ dfw1_k1; % starting time for Sw2_k1 at x = L1
        Xsw2_k1 = L1 + qt / A / phi2 .* dfw2_k1 .* (t2(ti) - ts); % travelling distance of Sw2_k1
        % 6.1.3 calculate saturation profiles for Swc2 <= Sw,2 <= Sw2,inf(1)
        krw2_k2 = kr_w2(Sw2_k2); % relative permeability of water
        kro2_k2 = kr_o2(Sw2_k2); % relative permeability of oil
        fw2_k2 = f_w2(Sw2_k2); % fractional flow function
        dfw2_k2 = df_w2(Sw2_k2); % derivate of fractional flow function
        Xsw2_k2 = L1 + qt / A / phi2 .* dfw2_k2 * (t2(ti) - t_inf); % travelling distance of Sw2_k2
        % 6.1.4 calculate the shock front saturation and its distance
        Xsw2_k{ti, 1} = [L1; Xsw2_k1; Xsw2_k2]; % the travelling distance of Sw2_k
        % delta_Xsw2 = Xsw2_j - Xsw2_(j-1), x2_0 = L1
        dlt_xsw2 = Xsw2_k{ti, 1}(2:end) - Xsw2_k{ti, 1}(1:end - 1);
        dlt_swk2 = Sw2_k{ti, 1}(2:end) - Swc2; % delta_Swk2 = Sw2,k - Swc2

        for i = 1:nsw2
            % calculate the injected fluid volume from L1 to Xswf
            V2 = A * phi2 * sum(dlt_xsw2(1:i) .* dlt_swk2(1:i));

            if V2 >= W2(ti)
                index2(ti) = i;
                Swf2(ti) = Sw2_k{ti, 1}(i);
                Xswf2(ti) = Xsw2_k{ti, 1}(i);
                Sw2_k{ti, 1}(i + 1:end) = Swc2;
                break;
            end

        end

        Xsw2_k{ti, 1} = Xsw2_k{ti, 1}(1:index2(ti));
        Sw2_k{ti, 1} = Sw2_k{ti, 1}(1:index2(ti));
    end

    % 6.2 figure 5-1: Plot the water saturation profiles in domain 1
    h_fig5 = figure(5);
    set(h_fig5, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Water Saturation Profiles');
    SwTime = [];
    subplot(1, 2, 1);

    for ti = 1:nts2
        plot(Xsw1(:, ti + nts1), Swt1, '-r', 'LineWidth', 2.0);
        hold on;
        plot([0 Xsw1(1, ti + nts1)], [1 - Sor1 1 - Sor1], '-r', 'LineWidth', 2.0);
        SwTime = [SwTime; ['Time = ' num2str(t2(ti) / 86400, format) ' day']];
    end

    plot([0 L1], [Sw2_inf(1) Sw2_inf(1)], '--k', 'LineWidth', 2.0);
    set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L1]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 14);
    title('domain 1');
    xlabel('{\it x} (m)'); ylabel('{\it S}_w');
    axis square;
    % 6.3 figure 5-2: Plot the water saturation profiles in domain 2
    subplot(1, 2, 2);

    for ti = 1:nts2
        plot(Xsw2_k{ti, 1}, Sw2_k{ti, 1}, '-b', 'LineWidth', 2.0);
        hold on;
        plot([Xsw2_k{ti, 1}(end) Xsw2_k{ti, 1}(end)], [Sw2_k{ti, 1}(end) Swc2], '-b',
        'LineWidth', 2.0);
        plot([Xsw2_k{ti, 1}(end) L1 + L2], [Swc2 Swc2], '-b', 'LineWidth', 2.0);
    end

    plot([L1 L1 + L2], [Sw2_inf(1) Sw2_inf(1)], '--k', 'LineWidth', 2.0);
    plot(L1, Sw2_inf(1), 'ro', 'markersize', 8);
    text(L1, Sw2_inf(1) - 0.05, '{\it S}_{w2}^{*}', 'Color', 'k', 'Fontname', 'Times New Roman', 'FontSize', 14);
    set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [L1 L1 + L2]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 14);
    set(gca, 'YAxisLocation', 'right');
    title('domain 2');
    xlabel('{\it x} (m)'); ylabel('{\it S}_w');
    h_legend5 = legend(SwTime);
    set(h_legend5, 'Box', 'on', 'Location', 'best');
    axis square;
    %% end

    % B.4 MATLAB CODE FOR BUCKLEY–LEVERETT SOLUTION
    % IN A RADIAL COMPOSITE SYSTEM
    %% Problem description
    % Analytical solution of Buckley–Leverett equation using Welge's method
    % in a radial composite porous medium system
    % Brooks–Corey type relative permeability curves are used
    close all; clear all; clc;

    %% 1 Parameters initialization for two domains
    % 1.1 rock properties in domain 1
    R1 = 4.00; % radius of domain 1 [m]
    phi1 = 0.30; % porosity of domain 1
    k1 = 1.0e-13; % permeability of domain 1 [m^2]
    Sor1 = 0.20; % residual oil saturation in domain 1
    Swc1 = 0.20; % connate water saturation in domain 1
    no1 = 1.50; % exponent of oil phase in domain
    nw1 = 2.50; % exponent of water phase in domain 1
    kro1_max = 0.75; % maximum permeability of oil domain 1
    krw1_max = 0.75; % maximum permeability of water domain 1
    dlt_Sw1 = 1e-3; % constant saturation step in domain 1
    Sw1 = [(Swc1 + eps):dlt_Sw1:(1 - Sor1 - eps)]'; % water saturation vector in domain 1
    % 1.2 rock properties in domain 2
    R2 = 8.00; % radius of domain 2 [m]
    phi2 = 0.30; % porosity of domain 2
    k2 = 1.0e-13; % permeability of domain 2 [m^2]
    Sor2 = 0.20; % residual oil saturation in domain 2
    Swc2 = 0.20; % connate water saturation in domain 2
    no2 = 2.50; % exponent of oil phase in domain 2
    nw2 = 1.50; % exponent of water phase in domain 2
    kro2_max = 0.80; % maximum permeability of oil domain 2
    krw2_max = 0.80; % maximum permeability of water domain 2
    dlt_Sw2 = 1e-3; % constant saturation step in domain 2
    Sw2 = [(Swc2 + eps):dlt_Sw2:(1 - Sor2 - eps)]'; % water saturation vector in domain 2
    % 1.3 fluid properties for the composite system
    muo = 5.0e-3; % oil phase viscosity [Pa*s]
    muw = 1.0e-3; % water phase viscosity [Pa*s]
    rho_o = 0.8e3; % oil density [kg/m^3]
    rho_w = 1.0e3; % water density [kg/m^3]
    dlt_rho = rho_w - rho_o; % density difference [kg/m^3]
    % 1.4 other parameters for calculation
    rw = 0.10; % radius of injection well [m]
    h = 1.0; % reservoir thickness [m]
    qt = 2.5e-4; % constant water injection rate [m^3/s]
    g = 9.8067; % gravity acceleration constant [m/s^2]
    time0 = 86400 * 0.1; % initial calculation time [s]
    timef = 86400 * 0.463; % final calculation time [s]
    nts1 = 5; % time steps for domain 1
    nts2 = 5; % time steps for domain 2
    nsw1 = size(Sw1, 1); % number of water saturation in domain 1
    nsw2 = size(Sw2, 1); % number of water saturation in domain 2
    format = '%4.3f'; % precision format for plotting

    %% 2 Calculate fractional flow function and plots
    % 2.1 initialization of function handles in domain 1
    kr_w1 = @(S) krw1_max .* ((S - Swc1) / (1 - Swc1 - Sor1)) .^ nw1; % relative permeability water
    kr_o1 = @(S) kro1_max .* ((1 - S - Sor1) / (1 - Swc1 - Sor1)) .^ no1; % relative permeability oil
    mob1 = @(S) (kr_o1(S) ./ muo) .* (muw ./ kr_w1(S)); % mobility function
    f_w1 = @(S) 1 ./ (1 + mob1(S)); % fractional flow function
    df_w1 = @(S) ((nw1 .* mob1(S) ./ (S - Swc1)) + (no1 .* mob1(S)) ./ (1 - S - Sor1))... % derivate
    ./ (1 + mob1(S)) .^ 2;
    % 2.2 initialization of function handles in domain 2
    kr_w2 = @(S) krw2_max .* ((S - Swc2) / (1 - Swc2 - Sor2)) .^ nw2; % relative permeability water
    kr_o2 = @(S) kro2_max .* ((1 - S - Sor2) / (1 - Swc2 - Sor2)) .^ no2; % relative permeability oil
    mob2 = @(S) (kr_o2(S) ./ muo) .* (muw ./ kr_w2(S)); % mobility function
    f_w2 = @(S) 1 ./ (1 + mob2(S)); % fractional flow function
    df_w2 = @(S) ((nw2 .* mob2(S) ./ (S - Swc2)) + (no2 .* mob2(S)) ./ (1 - S - Sor2))... % derivate
    ./ (1 + mob2(S)) .^ 2;
    % 2.3 relative permeability and fractional flow vectors in domain 1
    krw1 = kr_w1(Sw1); % relative permeability vector of water
    kro1 = kr_o1(Sw1); % relative permeability vector of oil
    fw1 = f_w1(Sw1); % fractional flow vector
    dfw1 = df_w1(Sw1); % fractional flow derivate vector
    % 2.4 relative permeability and fractional flow vectors in domain 2
    krw2 = kr_w2(Sw2); % relative permeability vector of water
    kro2 = kr_o2(Sw2); % relative permeability vector of oil
    fw2 = f_w2(Sw2); % fractional flow vector
    dfw2 = df_w2(Sw2); % fractional flow derivate vector
    % 2.5 plot relative permeability curves
    h_fig1 = figure(1);
    set(h_fig1, 'color', 'w', 'NumberTitle', 'off','Name', 'Relative Permeability Curves: kr');
    plot(Sw1, krw1, '-b', Sw1, kro1, '-r', 'LineWidth', 2.0);
    hold on;
    plot(Sw2, krw2, '--b', Sw2, kro2, '--r', 'LineWidth', 2.0);
    axis([0.0 1.0 0.0 1.0]);
    axis square;
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('Brooks–Corey relative permeability curves')
    xlabel('{\it S}_w');
    ylabel('{\it k}_{r{\it \beta}}');
    set(gca, 'YTick', 0:0.2:1);
    h_legend1 = legend('domain 1: water', 'domain 1: oil', 'domain 2: water', 'domain 2: oil');
    set(h_legend1, 'Box', 'on', 'Location', 'best');
    %% 3 Calculate advance frontal water saturation in domain 1 and plot
    % 3.1 calculate sock front of water saturation
    [index, Swf1, dfwf1, bf1] = calSatFront(Sw1, fw1, dfw1);
    % 3.2 Plot the fractional flow function curves
    h_fig2 = figure(2);
    set(h_fig2, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional Flow Curves: fw');
    plot(Sw1, fw1, '-b', 'LineWidth', 2.0);
    hold on;
    plot(Sw2, fw2, '--b', 'LineWidth', 2.0);
    plot(Sw1, (dfwf1 .* Sw1 + bf1), '-r', 'LineWidth', 2.0);
    plot(Swc1, 0, 'ro');
    plot(Swf1, fw1(index), 'ro');
    plot([Swf1 Swf1], [0 fw1(index)], '--r', 'LineWidth', 2.0);
    plot([0 Swf1], [fw1(index) fw1(index)], '--r', 'LineWidth', 2.0);
    axis([0.0 1.0 0.0 1.0]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('fractional flow function curves');
    xlabel('{\it S}_w');
    ylabel('{\it f}_w');
    h_legend2 = legend('domain 1', 'domain 2');
    set(gca, 'YTick', 0:0.2:1);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    set(h_legend2, 'Box', 'on', 'Location', 'best');
    axis square;
    % 3.3 Plot the derivative curve of fractional flow function
    h_fig3 = figure(3);
    set(h_fig3, 'color', 'w', 'NumberTitle', 'off','Name', 'Derivatives of Fractional Flow Curves: dfw/dSw');
    plot(Sw1, dfw1, '-b', 'LineWidth', 2.0);
    hold on;
    plot(Sw2, dfw2, '--b', 'LineWidth', 2.0);
    set(gca, 'XLim', [0 1]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('derivatives of fractional flow function');
    xlabel('{\it S}_w');
    ylabel('d {\it f}_w / d {\it S}_w');
    h_legend3 = legend('domain 1', 'domain 2');
    set(h_legend3, 'Box', 'on', 'Location', 'best');
    axis square;
    %% 4 Water saturation profiles in domain 1 before Swf1 reaches interface
    % 4.1 the time when Swf reaches at the interface between domain1 and domain 2
    krw1_swf = kr_w1(Swf1);
    kro1_swf = kr_o1(Swf1);
    fw1_swf = f_w1(Swf1);
    dfw1_swf = df_w1(Swf1);
    t_inf = pi * h * phi1 * (R1 ^ 2 - rw ^ 2) / (qt * dfw1_swf);
    % 4.2 calculate the travelling distance of specific water saturations Sw1
    t1 = linspace(time0, t_inf, nts1); % time vector for domain 1
    dfwt1 = dfw1(end:-1:index); % inverted sequence of dfw
    Swt1 = Sw1(end:-1:index); % inverted sequence of Sw

    for ti = 1:nts1 % loop times vector

        for i = 1:(nsw1 - index + 1)
            Rsw1(i, ti) = sqrt(rw ^ 2 + qt * t1(ti) / (pi * h * phi1) * dfwt1(i));
        end

    end

    % 4.3 figure 4: Plot the water saturation profiles in domain 1
    h_fig4 = figure(4);
    set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Water Saturation Profiles');
    SwTime = [];

    for ti = 1:nts1
        plot(Rsw1(:, ti), Swt1, '-r', 'LineWidth', 2.0);
        hold on;
        plot([rw Rsw1(1, ti)], [1 - Sor1 1 - Sor1], '-r', 'LineWidth', 2.0);
        SwTime = [SwTime; ['Time = ' num2str(t1(ti) / 86400, format) ' day']];
        plot([Rsw1(end, ti) Rsw1(end, ti)], [Swf1 Swc1], '-r', 'LineWidth', 2.0);
        plot([Rsw1(end, ti) R1], [Swc1 Swc1], '-r', 'LineWidth', 2.0);
    end

    set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 R1]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
    title('water saturation profiles in domain 1');
    xlabel('{\it x} (m)'); ylabel('{\it S}_w');
    h_legend4 = legend(SwTime);
    set(h_legend4, 'Box', 'on', 'Location', 'best');
    axis square;
    %% 5 Calculate water saturation profile in domain 1 when time > t_inf
    % 5.1 calculate the corresponding parameters at interface
    if timef <= t_inf
        error('Please reset the timef due to timef <= t_inf!');
    end

    t2 = linspace(t_inf, timef, nts2)';
    dfw1_inf = pi * h * phi1 * (R1 ^ 2 - rw ^ 2) / qt ./ t2; % dfw,1 at the interface (upstream)
    Sw1_inf = calObjFun(dfw1_inf, Sw1, dfw1, index); % Sw,1 at the interface (upstream)
    krw1_inf = kr_w1(Sw1_inf); % relative water permeability
    kro1_inf = kr_o1(Sw1_inf); % relative permeability of oil
    fw1_inf = f_w1(Sw1_inf); % fractional flow function for Sw1_inf
    Sw1_avg = Sw1_inf + (1 - fw1_inf) ./ dfw1_inf; % average saturation of domain 1
    % 5.2 % calculate the travelling distance of water saturations
    dfwt1 = dfw1(end:-1:index); % inverted sequence of dfw
    Swt1 = Sw1(end:-1:index); % inverted sequence of Sw

    for ti = 1:nts2 % loop times vector

        for i = 1:(nsw1 - index + 1)
            Rsw1(i, ti + nts1) = sqrt(rw ^ 2 + qt * t2(ti) / (pi * h * phi1) * dfwt1(i));
        end

    end

    %% 6 Calculate water saturation profile in domain 2 when time > t_inf
    % 6.1 calculate water saturation at interface and injected fluid
    fw2_inf = fw1_inf; % fw,2 at the interface (downstream)
    Sw2_inf = calObjFun(fw2_inf, Sw2, fw2); % Sw,2 at the interface (downstream)
    Qt = qt .* t2; % the injected water into domains
    W1 = pi * h * phi1 * (R1 ^ 2 - rw ^ 2) .* (Sw1_avg - Swc1); % the injected water into domain 1
    W2 = Qt - W1; % the injected water into domain 2
    if (abs(W2(1) / Qt(1)) <= 1e-3); W2(1) = 0; end; % enforce W2 == 0.0
    % 6.1 calculate saturation profiles using mass balance principle
    Sw2_k = cell(nts2, 1); % water saturation Sw,2: initialization
    Rsw2_k = cell(nts2, 1); % initial travelling distance of Sw,2
    Sw2_k{1, 1} = [Sw2_inf(1); Swc2]; % the saturation Sw,2 set at t_inf
    Rsw2_k{1, 1} = [R1; R1]; % the distance Xsw,2 set at t_inf

    for ti = 2:nts2
        % 6.1.1 select a set of saturation Sw,2
        % choose a set of saturation Sw,2: [Sw2,inf(1) Sw2,inf(ti)]
        Sw2_k1 = [(Sw2_inf(ti) - dlt_Sw2):(-dlt_Sw2):Sw2_inf(1)]';
        % choose a set of saturation Sw,2: [Swc2 Sw2,inf(1)])
        Sw2_k2 = [(Sw2_inf(1) - dlt_Sw2):(-dlt_Sw2):Swc2]';
        % the whole set of water saturation Sw,2
        Sw2_k{ti, 1} = [Sw2_inf(ti); Sw2_k1; Sw2_k2];
        % 6.1.2 calculate saturation profiles for Sw2,inf(1) <= Sw,2 <= Sw2,inf(ti)
        krw2_k1 = kr_w2(Sw2_k1); % relative permeability of water
        kro2_k1 = kr_o2(Sw2_k1); % relative permeability of oil
        fw2_k1 = f_w2(Sw2_k1); % fractional flow function
        dfw2_k1 = df_w2(Sw2_k1); % derivative of fractional flow function
        fw1_k1 = fw2_k1;
        Sw1_k1 = calObjFun(fw1_k1, Sw1, fw1); % calculate the Sw1 at interface
        krw1_k1 = kr_w1(Sw1_k1); % relative permeability of water
        kro1_k1 = kr_o1(Sw1_k1); % relative permeability of oil
        fw1_k1 = f_w1(Sw1_k1); % fractional flow function
        dfw1_k1 = df_w1(Sw1_k1); % derivative of fractional flow function
        ts = pi * h * phi1 * (R1 ^ 2 - rw ^ 2) / qt ./ dfw1_k1; % starting time for Sw,2 at r = R1
        % calculate the travelling distance of Sw2_k1
        Rsw2_k1 = sqrt(R1 ^ 2 + qt / (pi * h * phi2) .* dfw2_k1 .* (t2(ti) - ts));
        % 6.1.3 calculate saturation profiles for Swc2 <= Sw,2 <= Sw2,inf(1)
        krw2_k2 = kr_w2(Sw2_k2); % relative permeability of water
        kro2_k2 = kr_o2(Sw2_k2); % relative permeability of oil
        fw2_k2 = f_w2(Sw2_k2); % fractional flow function
        dfw2_k2 = df_w2(Sw2_k2); % derivative of fractional flow function
        % calculate the travelling distance of Sw2_k2
        Rsw2_k2 = sqrt(R1 ^ 2 + qt / (pi * h * phi2) .* dfw2_k2 .* (t2(ti) - t_inf));
        % 6.1.4 calculate the shock front saturation and its distance
        % calculate the travelling distance of Sw2_k
        Rsw2_k{ti, 1} = [R1; Rsw2_k1; Rsw2_k2];
        % delta_Rsw2 = Rsw2_j^2 - Rsw2_(j-1)^2, r2_0 = R1
        dlt_Rsw2 = Rsw2_k{ti, 1}(2:end) .^ 2 - Rsw2_k{ti, 1}(1:end - 1) .^ 2;
        % delta_Swk2 = Sw2,k - Swc2
        dlt_swk2 = Sw2_k{ti, 1}(2:end) - Swc2;

        for i = 1:nsw2
            % calculate the injected fluid volume from R1 to Xswf
            V2 = pi * h * phi2 * sum(dlt_Rsw2(1:i) .* dlt_swk2(1:i));

            if V2 >= W2(ti)
                index2(ti) = i;
                Swf2(ti) = Sw2_k{ti, 1}(i);
                Rswf2(ti) = Rsw2_k{ti, 1}(i);
                Sw2_k{ti, 1}(i + 1:end) = Swc2;
                break;
            end

        end

        Rsw2_k{ti, 1} = Rsw2_k{ti, 1}(1:index2(ti));
        Sw2_k{ti, 1} = Sw2_k{ti, 1}(1:index2(ti));
    end

    % 6.2 figure 5-1: Plot the water saturation profiles in domain 1
    h_fig5 = figure(5);
    set(h_fig5, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Water Saturation Profiles');
    SwTime = [];
    subplot(1, 2, 1);

    for ti = 1:nts2
        plot(Rsw1(:, ti + nts1), Swt1, '-r', 'LineWidth', 2.0);
        hold on;
        plot([rw Rsw1(1, ti + nts1)], [1 - Sor1 1 - Sor1], '-r', 'LineWidth', 2.0);
        SwTime = [SwTime; ['Time = ' num2str(t2(ti) / 86400, format) ' day']];
    end

    plot([rw R1], [Sw2_inf(1) Sw2_inf(1)], '--k', 'LineWidth', 2.0);
    set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 R1]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 14);
    title('domain 1');
    xlabel('{\it x} (m)'); ylabel('{\it S}_w');
    axis square;
    % 6.3 figure 5-2: Plot the water saturation profiles in domain 2
    subplot(1, 2, 2);

    for ti = 1:nts2
        plot(Rsw2_k{ti, 1}, Sw2_k{ti, 1}, '-b', 'LineWidth', 2.0);
        hold on;
        plot([Rsw2_k{ti, 1}(end) Rsw2_k{ti, 1}(end)], [Sw2_k{ti, 1}(end) Swc2], '-b','LineWidth', 2.0);
        plot([Rsw2_k{ti, 1}(end) R1 + R2], [Swc2 Swc2], '-b', 'LineWidth', 2.0);
    end

    plot([R1 R2], [Sw2_inf(1) Sw2_inf(1)], '–k', 'LineWidth', 2.0);
    plot(R1, Sw2_inf(1), 'ro', 'markersize', 8);
    text(R1, Sw2_inf(1) - 0.05, '{\it S}_{w2}^{*}', 'Color', 'k', 'Fontname','Times New Roman', 'FontSize', 14);
    set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [R1 R2]);
    set(gca, 'Fontname', 'Times New Roman', 'FontSize', 14);
    set(gca, 'YAxisLocation', 'right');
    title('domain 2');
    xlabel('{\it x} (m)'); ylabel('{\it S}_w');
    h_legend5 = legend(SwTime);
    set(h_legend5, 'Box', 'on', 'Location', 'best');
    axis square;
    %% end

    %B.5 AUXILIARY FUNCTIONS
    function [index, Swf, dfwf, bf] = calSatFront(Sw, fw, dfw)
        %% calculate advance front saturation using Welge graphic method
        % Sw ––– saturation vector of displacing wetting phase fluid
        % fw ––– fractional flow vector of displacing wetting phase fluid
        % dfw ––– fractional flow derivate vector of displacing wetting phase fluid
        % index ––– index of Swf in Sw
        % Swf ––– advance front saturation
        % dfwf ––– fractional flow derivative at Swf
        % bf ––– intercept valuve for tangent line
        index = 0;
        ns = size(Sw, 1);
        df = 1 / (Sw(ns) - Sw(1));
        dlt_Sw = Sw(ns) - Sw(ns - 1);

        if nargin == 2
            % calculate derivate
            for i = 3:(ns - 2)
                dfw(i) = (fw(i + 1) - fw(i - 1)) / (2 * dlt_Sw);
            end

            dfw(2) = (-11 * fw(2) + 18 * fw(3) - 9 * fw(4) + 2 * fw(5)) / (6 * dlt_Sw);
            dfw(1) = abs(2 * dfw(2) - dfw(3));
            dfw(ns - 1) =- (-11 * fw(nt - 1) + 18 * fw(nt - 2) - 9 * fw(nt - 3) + 2 * fw(nt - 4)) / (6 * dlt_Sw);
            dfw(ns) = abs(2 * dfw(nt - 1) - dfw(nt - 2));
            % calculate Swf and its index
            for i = 2:ns - 1
                dfds = fw(i) / (Sw(i) - Sw(1));

                if (dfds < dfw(i - 1)) && (dfds > dfw(i + 1)) && (dfds >= df)
                    index = i;
                    Swf = Sw(i);
                    dfwf = dfds;
                    bf = fw(i) - dfwf * Sw(i);
                    break;
                end

            end

            if index == 0
                error('Can''t find a tagent line through point (Swc,0). Decrease dlt_Sw!');
            end

        elseif nargin == 3
            % calculate Swf and its index
            for i = 2:ns - 1
                dfds = fw(i) / (Sw(i) - Sw(1));

                if (dfds < dfw(i - 1)) && (dfds > dfw(i + 1)) && (dfds >= df)
                    index = i;
                    Swf = Sw(i);
                    dfwf = dfds;
                    bf = fw(i) - dfwf * Sw(i);
                    break;
                end

            end

            if index == 0
                error('Can''t find a tagent line through point (Swc,0). Decrease dlt_Sw!');
            end

        else
            error('the number of input arguments in function calSatFront is incorrect!');
        end

    end

    %% end function calSatFront()
    function obj_fun = calObjFun(obj_var, Fun, Var, index_var)
        %% calculate obj_fun at obj_var using interpolation method
        % Fun ––– basic function value vector
        % Var ––– basic variable value vector
        % obj_var ––– objective variable value vector
        % obj_fun ––– objective function value vector
        % index ––– optional for non-monotone Fun
        ns = size(obj_var, 1);
        obj_fun = zeros(ns, 1);

        if nargin == 3

            for i = 1:ns
                [w, index] = sort(abs(Var - obj_var(i)));
                obj_fun(i) = Fun(index(1));
            end

        elseif nargin == 4

            for i = 1:ns
                [w, index] = sort(abs(Var - obj_var(i)));

                if index(1) >= index_var
                    obj_fun(i) = Fun(index(1));
                else
                    obj_fun(i) = Fun(index(2));
                end

            end

        else
            error('Wrong input parameters in calObjFun function!');
        end

    end

    %% end function calObjFun()
