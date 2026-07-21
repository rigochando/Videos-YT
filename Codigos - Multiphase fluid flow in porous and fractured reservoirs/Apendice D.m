%D.1 MATLAB CODE FOR BUCKLEY–LEVERETT SOLUTION:
%FORCHHEIMER EQUATION

%% Problem description
% Analytical solution of Buckley–Leverett equation using integral method
% non-Darcy flow: wetting phase fluid displacing non-wetting phase fluid
% non-Darcy flow is described by Forchheimer equation
% Brooks-Corey relative permeability curve has been applied
close all; clear all; clc;
%% 1 Parameters initialization
% 1.1 rock properties
L = 5.00; % domain length [m]
phi = 0.30; % porosity
k = 9.869e-12; % absolute permeability [m^2]
A = 1.0; % area of cross-section [m^2]
theta = pi * 0.0; % angle between x and horizontal [rad]
% 1.2 fluid properties
mun = 5.0e-3; % non-wetting phase viscosity [Pa*s]
muw = 1.0e-3; % wetting phase viscosity [Pa*s]
rho_n = 0.8e3; % non-wetting phase density [kg/m^3]
rho_w = 1.0e3; % wetting phase density [kg/m^3]
dlt_rho = rho_w - rho_n; % phase density difference [kg/m^3]
% 1.3 relative permeability parameters
Snr = 0.20; % residual saturation non-wetting phase
Swr = 0.20; % residual saturation wetting phase
nn = 2.00; % exponent of non-wetting phase
nw = 2.00; % exponent of wetting phase
krn_max = 0.75; % maximum permeability of non-wetting
krw_max = 0.75; % maximum permeability of wetting phase
% 1.4 other initial parameters
C_bw = 3.20e-6; % non-Darcy flow constant [m^(3/2)]
C_bn = 3.20e-6; % non-Darcy flow constant [m^(3/2)]
dlt_Sw = 1e-3; % constant saturation step
Sw = [(Swr+eps):dlt_Sw:(1 - Snr - eps)]'; % wetting phase saturation vector
qt = 5.0e-4; % injection rate [m^3/s]
px0 = -1e5; % initial test pressure gradient [Pa/m]
g = 9.80665; % gravity acceleration constant [m/s^2]
time0 = 3600 * 0.02; % initial calculation time [s]
timef = 3600 * 0.20; % final calculation time [s]
nts = 4; % time steps for calculating
nsw = size(Sw, 1); % size of Sw vector
format = '%3.2e'; % precision format for plotting
%% 2 Relative permeabilities, fractional flow function and its derivative
% 2.1 initialization of function handles
kr_w = @(S) krw_max .* ((S - Swr) / (1 - Swr - Snr)) .^ nw; % relative permeability function
kr_n = @(S) krn_max .* ((1 - S - Snr) / (1 - Swr - Snr)) .^ nn; % relative permeability function
b_w = @(S) C_bw ./ (k .* kr_w(S)) .^ (5/4) ./ (phi .* (S - Swr)) .^ (3/4); % coefficient beta_w
b_n = @(S) C_bn ./ (k .* kr_n(S)) .^ (5/4) ./ (phi .* (1 - S - Snr)) .^ (3/4); % coefficient beta_n
v_w = @(S, p_x) 1 ./ (2 * k * rho_w .* b_w(S)) .* (-muw ./ kr_w(S)+ ... % wetting phase
    sqrt((muw ./ kr_w(S)) .^ 2 - 4 * k ^ 2 * rho_w .* b_w(S) .* (p_x+rho_w * g * sin(theta))));
v_n = @(S, p_x) 1 ./ (2 * k * rho_n .* b_n(S)) .* (-mun ./ kr_n(S)+ ... % non-wetting phase
    sqrt((mun ./ kr_n(S)) .^ 2 - 4 * k ^ 2 * rho_n .* b_n(S) .* (p_x+rho_n * g * sin(theta))))
% 2.2 calculate the pressure gradient at specific saturation
px = zeros(nsw, 1); % initialization of dp/dx

for i = 1:nsw
    F = @(p_x) qt / A - v_w(Sw(i), p_x) - v_n(Sw(i), p_x);
    px(i) = fzero(F, px0);

    if isnan(px(i))
        error('Can not find the zero near x0! Reset parameters.');
    end

    px0 = px(i);
end

% 2.3 calculate the velocity, fractional flow function and its derivate
krw = kr_w(Sw); % relative permeability vector
krn = kr_n(Sw); % relative permeability vector
bw = b_w(Sw); % non-Darcy flow coefficient beta_w
bn = b_n(Sw); % non-Darcy flow coefficient beta_n
vw = v_w(Sw, px); % wetting phase velocity [m/s]
vn = v_n(Sw, px); % non-wetting phase velocity [m/s]
fw = vw ./ (vw+vn); % fractional flow function
dfw = zeros(nsw, 1); % derivate of fractional flow function

for i = 3:(nsw - 2)
    dfw(i) = (fw(i+1) - fw(i - 1)) / (2 * dlt_Sw);
end

dfw(2) = (-11 * fw(2)+18 * fw(3) - 9 * fw(4)+2 * fw(5)) / (6 * dlt_Sw);
dfw(1) = abs(2 * dfw(2) - dfw(3));
dfw(nsw - 1) =- (-11 * fw(nsw - 1)+18 * fw(nsw - 2) - 9 * fw(nsw - 3)+2 * fw(nsw - 4)) / (6 * dlt_Sw);
dfw(nsw) = abs(2 * dfw(nsw - 1) - dfw(nsw - 2));
%% 3 Calculate advance frontal saturation and travelling distance
% 3.1 calculate travelling distance of every saturations
t = linspace(time0, timef, nts);
dfwt = dfw(end:-1:1); % inverted sequence of dfw
fwt = fw(end:-1:1); % inverted sequence of fw
Swt = Sw(end:-1:1); % inverted sequence of Sw

for ti = 1:nts
    Wi(ti) = qt * t(ti); % injected wetting phase fluid volume

    for i = 1:nsw
        Xsw(i, ti) = qt * t(ti) / (A * phi) * dfwt(i); % calculate travelling distance
    end

    Xsw0 = [0; Xsw(1:end - 1, ti)]; % the 2nd travelling distance vector
    dlt_Xsw = Xsw(:, ti) - Xsw0; % dlt_Xsw = Xsw_j - Xsw_(j-1), x_0 = 0
    dlt_Swt = Swt - Swr; % dlt_Swt = Swt,j - Swr

    for i = 1:nsw
        % calculate the injected fluid volume from 0 to Xswf
        Vi = A * phi * sum(dlt_Xsw(1:i) .* dlt_Swt(1:i));

        if Vi >= Wi(ti)
            index(ti) = i; % index of the Swf in Swt
            Swf(ti) = Swt(i); % advance front saturation: Swf
            Xswf(ti) = Xsw(i, ti); % the position of Swf
            break;
        end

    end

end

Xsw = Xsw(1:index(1), :);
Swt = Swt(1:index(1));
% 3.2 Calculate the time when Swf reaches at production well
krw_swf = kr_w(Swf(1));
krn_swf = kr_n(Swf(1));
fw_swf = calObjFun(Swf(1), fw, Sw);
dfw_swf = calObjFun(Swf(1), dfw, Sw);
t_pw = A * phi * L / (qt * dfw_swf);

%% 4 Plot results
% 4.1 Plot the relative permeability curves
h_fig1 = figure(1);
set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves:kr');
plot(Sw, krw, '-b', Sw, krn, '-r', 'LineWidth', 2.0);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Brooks-Corey relative permeability curves')
xlabel('{\it S}_w');
ylabel('{\it k}_{r {\it l}}');
set(gca, 'YTick', 0:0.2:1);
h_legend1 = legend('wetting phase', 'non-wetting phase');
set(h_legend1, 'Box', 'on', 'Location', 'best');

% 4.2 figure 2: Plot the dp/dx–Sne curves
px = -px / 1e5;
h_fig2 = figure(2);
set(h_fig2, 'color', 'w', 'NumberTitle', 'off', ...
    'Name', 'Pressure gradient versus Sw');
plot(Sw, px, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('pressure gradient');
xlabel('{\it S}_w');
ylabel('-d {\it p} / d {\it x}');
axis square;

% 4.3 Plot the non-Darcy flow coefficients curve
h_fig3 = figure(3);
set(h_fig3, 'color', 'w', 'NumberTitle', 'off', 'Name', 'non-Darcy coefficients');
semilogy(Sw, bw, '-b', Sw, bn, '-r', 'LineWidth', 2.0);
axis square;
axis([0.0 1.0 0.0 10 ^ 18]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('non-Darcy flow coefficient')
xlabel('{\it S}_w');
ylabel('{\it \betal}_{l} (m^{-1})');
h_legend3 = legend('wetting phase', 'non-wetting phase');
set(h_legend3, 'Box', 'on', 'Location', 'best');

% 4.4 Plot the fractional flow function curve
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional Flow Curve fw');
plot(Sw, fw, '-b', 'LineWidth', 2.0);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('fractional flow function curve')
xlabel('{\it S}_w');
ylabel('{\it f}_w');
set(gca, 'YTick', 0:0.2:1);

% 4.5 Plot the derivative curve of fractional flow function
h_fig5 = figure(5);
set(h_fig5, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Derivate of Fractional Flow Curve dfw / dSw');
plot(Sw, dfw, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('derivatives of fractional flow function')
xlabel('{\it S}_w');
ylabel('d {\it f}_w / d {\it S}_w');
axis square;
% 4.6 Plot the wetting phase saturation profiles
h_fig6 = figure(6);
set(h_fig6, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Saturation Profiles: Sw(t)');
SwTime = [];

for ti = 1:nts
    plot(Xsw(:, ti), Swt, '-r', 'LineWidth', 2.0);
    hold on;
    plot([0 Xsw(1, ti)], [1 - Snr 1 - Snr], '-r', 'LineWidth', 2.0);
    SwTime = [SwTime; ['Time = ' num2str(t(1, ti) / 3600, format) ' hours']];
    plot([Xsw(end, ti) Xsw(end, ti)], [Swf(ti) Swr], '-r', 'LineWidth', 2.0);
    plot([Xsw(end, ti) L], [Swr Swr], '-r', 'LineWidth', 2.0);
end

set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('wetting phase saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_w');
h_legend6 = legend(SwTime);
set(h_legend6, 'Box', 'on', 'Location', 'best');
axis square;
%% end

%D.2 MATLAB CODE FOR BUCKLEY–LEVERETT SOLUTION:
%BARREE AND CONWAY MODEL

%% Problem description
% Analytical solution of Buckley–Leverett equation using integral method
% non-Darcy flow: wetting phase fluid displacing non-wetting phase fluid
% non-Darcy flow is described by Barree–Conway model
% Brooks-Corey relative permeability curve has been applied
close all; clear all; clc;
%% 1 Parameters initialization
% 1.1 rock properties
L = 5.00; % domain length [m]
phi = 0.30; % porosity
kd = 9.869e-13; % absolute permeability [m^2]
kmr = 0.01; % minimum permeability ratio (fraction)
tau = 5e2; % inverse of characteristic length [1/m]
A = 1.0; % area of cross-section [m^2]
theta = pi * 0.0; % angle between x and horizontal [rad]
% 1.2 fluid properties
mun = 5.0e-3; % non-wetting phase viscosity [Pa*s]
muw = 1.0e-3; % wetting phase viscosity [Pa*s]
rho_n = 0.8e3; % non-wetting phase density [kg/m^3]
rho_w = 1.0e3; % wetting phase density [kg/m^3]
dlt_rho = rho_w - rho_n; % phase density difference [kg/m^3]
% 1.3 relative permeability parameters
Snr = 0.20; % residual saturation non-wetting phase
Swr = 0.20; % residual saturation wetting phase
nn = 2.00; % exponent of non-wetting phase
nw = 2.00; % exponent of wetting phase
krn_max = 0.75; % maximum permeability non-wetting phase
krw_max = 0.75; % maximum permeability wetting phase
% 1.4 other initial parameters
dlt_Sw = 1e-3; % constant saturation step
Sw = [(Swr+eps):dlt_Sw:(1 - Snr - eps)]'; % wetting phase saturation vector
px0 = -5e5; % find a zero of function near x_0
qt = 1.0e-5; % injection rate [m^3/s]
g = 9.80665; % gravity acceleration constant [m/s^2]
time0 = 3600 * 1.00; % initial calculation time [s]
timef = 3600 * 10.0; % final calculation time [s]
nts = 4; % time steps for calculating
nsw = size(Sw, 1); % size of Sw vector
format = '%3.2e'; % precision format for plotting
%% 2 Relative permeabilities, fractional flow function and its derivative
% 2.1 initialization of function handles
kr_w = @(S) krw_max .* ((S - Swr) / (1 - Swr - Snr)) .^ nw; % relative permeability function
kr_n = @(S) krn_max .* ((1 - S - Snr) / (1 - Swr - Snr)) .^ nn; % relative permeability function
v_w = @(S, p_x) -1 ./ (2 * muw * rho_w) .* ... % non-Darcy flow velocity
    (muw ^ 2 .* S * tau+(p_x+rho_w * g * sin(theta)) * kd .* kr_w(S) * kmr * rho_w)...
+1 ./ (2 * muw * rho_w) .* ...
    ((muw ^ 2 .* S * tau+(p_x+rho_w * g * sin(theta)) * kd .* kr_w(S) * kmr * rho_w) .^ 2 - ...
    4 * muw ^ 2 * rho_w * kd * tau .* S .* kr_w(S) .* (p_x+rho_w * g * sin(theta))) .^ (1/2);
v_n = @(S, p_x) -1 ./ (2 * mun * rho_n) .* ... % non-Darcy flow velocity
    (mun ^ 2 .* (1 - S) * tau+(p_x+rho_n * g * sin(theta)) * kd .* kr_n(S) * kmr * rho_n) +1 ./ (2 * mun * rho_n) .* ...
    ((mun ^ 2 .* (1 - S) * tau+(p_x+rho_n * g * sin(theta)) * kd .* kr_n(S) * kmr * rho_n) .^ 2 - ...
    4 * mun ^ 2 * rho_n * kd * tau .* (1 - S) .* kr_n(S) .* (p_x+rho_n * g * sin(theta))) .^ (1/2);
% 2.2 calculate the pressure gradient at specific saturation
px = zeros(nsw, 1); % initialization of dp/dx

for i = 1:nsw
    F = @(p_x) qt / A - v_w(Sw(i), p_x) - v_n(Sw(i), p_x);
    px(i) = fzero(F, px0);

    if isnan(px(i))
        error('Can not find the px near px0! Reset parameters.');
    end

    px0 = px(i);
end

% 2.3 calculate the velocity, fractional flow function and it's derivate
krw = kr_w(Sw); % relative permeability wetting phase
krn = kr_n(Sw); % relative permeability non-wetting phase
vw = v_w(Sw, px); % wetting phase velocity vector [m/s]
vn = v_n(Sw, px); % nonwetting phase velocity vector [m/s]
fw = vw ./ (vw+vn); % fractional flow function
dfw = zeros(nsw, 1); % derivate of fractional flow function

for i = 3:(nsw - 2)
    dfw(i) = (fw(i+1) - fw(i - 1)) / (2 * dlt_Sw);
end

dfw(2) = (-11 * fw(2)+18 * fw(3) - 9 * fw(4)+2 * fw(5)) / (6 * dlt_Sw);
dfw(1) = abs(2 * dfw(2) - dfw(3));
dfw(nsw - 1) =- (-11 * fw(nsw - 1)+18 * fw(nsw - 2) - 9 * fw(nsw - 3)+2 * fw(nsw - 4)) / (6 * dlt_Sw);
dfw(nsw) = abs(2 * dfw(nsw - 1) - dfw(nsw - 2));

%% 3 Calculate advance frontal saturation and travelling distance
% 3.1 calculate travelling distance of every saturations
t = linspace(time0, timef, nts);
dfwt = dfw(end:-1:1); % inverted sequence of dfw
fwt = fw(end:-1:1); % inverted sequence of fw
Swt = Sw(end:-1:1); % inverted sequence of Sw

for ti = 1:nts
    Wi(ti) = qt * t(ti); % injected wetting phase fluid volume

    for i = 1:nsw
        Xsw(i, ti) = qt * t(ti) / (A * phi) * dfwt(i); % calculate travelling distance
    end

    Xsw0 = [0; Xsw(1:end - 1, ti)]; % the 2nd travelling distance vector
    dlt_Xsw = Xsw(:, ti) - Xsw0; % dlt_Xsw = Xsw_j - Xsw_(j-1), x_0 = 0
    dlt_Swt = Swt - Swr; % dlt_Swt = Swt,j - Swr

    for i = 1:nsw
        % calculate the injected fluid volume from 0 to Xswf
        Vi = A * phi * sum(dlt_Xsw(1:i) .* dlt_Swt(1:i));

        if Vi >= Wi(ti)
            index(ti) = i; % index of the Swf in Swt
            Swf(ti) = Swt(i); % advance front saturation: Swf
            Xswf(ti) = Xsw(i, ti); % the position of Swf
            break;
        end

    end

end

Xsw = Xsw(1:index(1), :);
Swt = Swt(1:index(1));
% 3.2 Calculate the time when Swf reaches at production well
krw_swf = kr_w(Swf(1));
krn_swf = kr_n(Swf(1));
fw_swf = calObjFun(Swf(1), fw, Sw);
dfw_swf = calObjFun(Swf(1), dfw, Sw);
t_pw = A * phi * L / (qt * dfw_swf);
%% 4 Plot results
% 4.1 Plot the relative permeability curves
h_fig1 = figure(1);
set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves:kr');
plot(Sw, krw, '-b', Sw, krn, '-r', 'LineWidth', 2.0);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('(Brooks-Corey relative permeability curves')
xlabel('{\it S}_w');
ylabel('{\it k}_{r {\it l}}');
set(gca, 'YTick', 0:0.2:1);
h_legend1 = legend('wetting phase', 'non-wetting phase');
set(h_legend1, 'Box', 'on', 'Location', 'best');

% 4.2 figure 2: Plot the dp/dx––Sne curves
px = -px / 1e5;
h_fig2 = figure(2);
set(h_fig2, 'color', 'w', 'NumberTitle', 'off', ...
    'Name', 'pressure gradient curve');
plot(Sw, px, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('pressure gradient');
xlabel('{\it S}_{ne}');
ylabel('-d {\it p} / d {\it x}');
axis square;

% 4.3 Plot the fractional flow function curve
h_fig3 = figure(3);
set(h_fig3, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional Flow Curve fw');
plot(Sw, fw, '-b', 'LineWidth', 2.0);
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('fractional flow function curve')
xlabel('{\it S}_w');
ylabel('{\it f}_w');
set(gca, 'YTick', 0:0.2:1);

% 4.4 Plot the derivative curve of fractional flow function
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Derivates of Fractional Flow Curve dfw / dSw');
plot(Sw, dfw, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('derivate of fractional flow function')
xlabel('{\it S}_w');
ylabel('d {\it f}_w / d {\it S}_w');
axis square;
% 4.5 Plot the water saturation profiles
h_fig5 = figure(5);
set(h_fig5, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Saturation Profiles: Sw(t)');
SwTime = [];

for ti = 1:nts
    plot(Xsw(:, ti), Swt, '-r', 'LineWidth', 2.0);
    hold on;
    plot([0 Xsw(1, ti)], [1 - Snr 1 - Snr], '-r', 'LineWidth', 2.0);
    SwTime = [SwTime; ['Time = ' num2str(t(1, ti) / 3600, format) ' hours']];
    plot([Xsw(end, ti) Xsw(end, ti)], [Swf(ti) Swr], '-r', 'LineWidth', 2.0);
    plot([Xsw(end, ti) L], [Swr Swr], '-r', 'LineWidth', 2.0);
end

set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('wetting phase saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_w');
h_legend5 = legend(SwTime);
set(h_legend5, 'Box', 'on', 'Location', 'best');
axis square;
%% end

%D.3 AUXILIARY FUNCTIONS
function [index, Swf, dfwf, bf] = calSatFront(Sw, fw, dfw)
    %% calculate advance front saturation using Welge graphic method
    % Sw––– saturation vector of displacing wetting phase fluid
    % fw––– fractional flow vector of displacing wetting phase fluid
    % dfw ––– fractional flow derivate vector of displacing wetting phase fluid
    % index ––– index of Swf in Sw
    % Swf ––– advance front saturation
    % dfwf ––– fractional flow derivative at Swf
    % bf––– intercept value for tangent line
    index = 0;
    ns = size(Sw, 1);
    df = 1 / (Sw(ns) - Sw(1));
    dlt_Sw = Sw(ns) - Sw(ns - 1);

    if nargin == 2
        % calculate derivate
        for i = 3:(ns - 2)
            dfw(i) = (fw(i+1) - fw(i - 1)) / (2 * dlt_Sw);
        end

        dfw(2) = (-11 * fw(2)+18 * fw(3) - 9 * fw(4)+2 * fw(5)) / (6 * dlt_Sw);
        dfw(1) = abs(2 * dfw(2) - dfw(3));
        dfw(ns - 1) =- (-11 * fw(nt - 1)+18 * fw(nt - 2) - 9 * fw(nt - 3)+2 * fw(nt - 4)) / (6 * dlt_Sw);
        dfw(ns) = abs(2 * dfw(nt - 1) - dfw(nt - 2));
        % calculate Swf and its index
        for i = 2:ns - 1
            dfds = fw(i) / (Sw(i) - Sw(1));

            if (dfds < dfw(i - 1)) && (dfds > dfw(i+1)) && (dfds >= df)
                index = i;
                Swf = Sw(i);
                dfwf = dfds;
                bf = fw(i) - dfwf * Sw(i);
                break;
            end

        end

        if index == 0
            error('Cant find a tangent line through point (Swc, 0). Decrease dlt_Sw !');
        end

    elseif nargin == 3
        % calculate Swf and its index
        for i = 2:ns - 1
            dfds = fw(i) / (Sw(i) - Sw(1));

            if (dfds < dfw(i - 1)) && (dfds > dfw(i+1)) && (dfds >= df)
                index = i;
                Swf = Sw(i);
                dfwf = dfds;
                bf = fw(i) - dfwf * Sw(i);
                break;
            end

        end

        if index == 0
            error('Cant find a tangent line through point (Swc, 0). Decrease dlt_Sw !');
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
