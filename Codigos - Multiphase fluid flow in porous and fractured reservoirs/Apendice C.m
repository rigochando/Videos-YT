%C.1 MATLAB CODE FOR BUCKLEY–LEVERETT SOLUTION:
%POWER-LAW NON-NEWTONIAN FLUID

%% Problem description
% Analytical solution of Buckley–Leverett equation using Welge method
% A Newtonian fluid displaced by a Power-law non-Newtonian fluid
% Brooks-Corey relative permeability curve has been applied
close all; clear all; clc;
%% 1 Parameters initialization
% 1.1 rock properties
L = 5.00; % domain length [m]
phi = 0.20; % porosity
k = 9.869e-13; % absolute permeability [m^2]
A = 1.0; % area of cross-section [m^2]
alpha = pi * 0.0; % angle between x and horizontal [rad]
% 1.2 fluid properties
mune = 6.0e-3; % Newtonian phase viscosity [Pa*s]
rho_ne = 1.0e3; % Newtonian phase density [kg/m^3]
rho_nn = 0.8e3; % non-Newtonian phase density [kg/m^3]
dlt_rho = rho_nn - rho_ne; % phase density difference [kg/m^3]
n = 0.60; % power-law index
H = 0.010; % power-law coefficient [Pa*s^n]
% 1.3 relative permeability parameters
Sner = 0.15; % residual Newtonian saturation
Snnr = 0.15; % residual non-Newtonian saturation
nn = 3.00; % exponent of non-Newtonian phase
ne = 3.00; % exponent of Newtonian phase
krnn_max = 0.85; % maximum permeability of non-Newtonian
krne_max = 0.85; % maximum permeability of Newtonian
% 1.4 other initial parameters
dlt_Snn = 1e-4; % constant saturation step
eps = 1e-5; % a very small number
Snn = [(Snnr + eps):dlt_Snn:(1 - Sner - eps)]'; % wetting phase saturation vector
qt = 1.0e-5; % injection rate [m^3/s]
px0 = -5.5e4; % initial pressure gradient, px [Pa/m]
g = 9.80665; % gravity acceleration constant [m/s^2]
time0 = 3600 * 1.00; % initial calculation time [s]
timef = 3600 * 10.0; % final calculation time [s]
nts = 4; % time steps for calculating
nsn = size(Snn, 1); % size of Snn
format = '%3.2e'; % precision format for plotting
%% 2 Relative permeabilities, fractional flow function and its derivative
kr_nn = @(S) krnn_max .* ((S - Snnr) ./ (1 - Snnr - Sner)) .^ nn; % relative permeability
kr_ne = @(S) krne_max .* ((1 - S - Sner) ./ (1 - Snnr - Sner)) .^ ne; % relative permeability
mu_ef = @(S) H / 12 * (9 + 3 / n) ^ n * ... % effective viscosity
(150 * k .* kr_nn(S) * phi .* (S - Snnr)) .^ ((1 - n) / 2);
mu_nn = @(S, p_x) ... % apparent viscosity
mu_ef(S) .* (k .* kr_nn(S) ./ mu_ef(S) .* abs(p_x)) .^ (1 - 1 / n);
v_nn = @(S, p_x) -k .* kr_nn(S) / mu_nn(S, p_x) .* ... % velocity: nn
(p_x + rho_nn * g * sin(alpha));
v_ne = @(S, p_x) -k .* kr_ne(S) / mune .* (p_x + rho_ne * g * sin(alpha)); % velocity: ne
px = zeros(nsn, 1); % initialization of dp/dx

for i = 1:nsn
    F = @(p_x) qt / A - v_nn(Snn(i), p_x) - v_ne(Snn(i), p_x);
    px(i) = fzero(F, px0);

    if isnan(px(i))
        error('Can not find the zero near x0! Reset parameters.');
    end

    px0 = px(i);
end

krnn = kr_nn(Snn); % relative permeability: non-Newtonian
krne = kr_ne(Snn); % relative permeability: Newtonian fluid
muef = mu_ef(Snn); % effective viscosity: non-Newtonian
munn = mu_nn(Snn, px); % non-Newtonian fluid viscosity
mob = (krne ./ mune) .* (munn ./ krnn); % mobility function
fnn = 1 ./ (1 + mob) - ... % fractional flow function
A * k .* krne ./ (qt * mune) * dlt_rho * g * sin(alpha) ./ (1 + mob);
dfnn = zeros(nsn, 1); % derivate of fractional flow function

for i = 3:(nsn - 2)
    dfnn(i) = (fnn(i + 1) - fnn(i - 1)) / (2 * dlt_Snn);
end

dfnn(2) = (-11 * fnn(2) + 18 * fnn(3) - 9 * fnn(4) + 2 * fnn(5)) / (6 * dlt_Snn);
dfnn(1) = abs(2 * dfnn(2) - dfnn(3));
dfnn(nsn - 1) =- (-11 * fnn(nsn - 1) + 18 * fnn(nsn - 2) - 9 * fnn(nsn - 3) + 2 * fnn(nsn - 4)) /(6 * dlt_Snn);
dfnn(nsn) = abs(2 * dfnn(nsn - 1) - dfnn(nsn - 2));
% figure 1: Plot the relative permeability curves
h_fig1 = figure(1);
set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves: kr');
plot(Snn, krnn, '-b', Snn, krne, '-r', 'LineWidth', 2.0);
hold on;
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Brooks-Corey relative permeability curves')
xlabel('{\it S}_{nn}');
ylabel('{\it k}_{r{\it \beta}}');
set(gca, 'YTick', 0:0.2:1);
h_legend1 = legend('non-Newtonian phase', 'Newtonian phase');
set(h_legend1, 'Box', 'on', 'Location', 'best');
% figure 2: Plot the dp/dx--Snn curves
h_fig2 = figure(2);
set(h_fig2, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Pressure gradient versus Snn');
plot(Snn, -px / 1e5, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Pressure gradient (Bar/m)');
xlabel('{\it S}_{nn}');
ylabel('-d {\it p} / d {\it x}');
axis square;
% figure 2: Plot the fractional flow function curves
h_fig3 = figure(3);
set(h_fig3, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional flow function curve');
plot(Snn, fnn, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Fractional flow function');
xlabel('{\it S}_{nn}');
ylabel('{\it f}_{nn}');
axis square;

%% 3 Calculate advance frontal saturation
% 3.1 Evaluate advance frontal non-Newtonian fluid saturation Snnf
[index, Snnf, dfnnf, bf] = calSatFront(Snn, fnn, dfnn);
% 3.2 Evaluate the time of advance frontal saturation Snnf
krnn_sf = kr_nn(Snnf);
krne_sf = kr_ne(Snnf);
munn_sf = calObjFun(Snnf, munn, Snn);
mob_sf = (munn_sf / krnn_sf) * (krne_sf / mune); % mobility function
fnn_sf = 1 ./ (1 + mob_sf) - ... % fractional flow function
A * k .* krne_sf ./ (qt * mune) * dlt_rho * g * sin(alpha) ./ (1 + mob_sf);
dfnn_sf = calObjFun(Snnf, dfnn, Snn); % derivative of fractional flow function
tpw = A * phi * L / (qt * dfnn_sf); % time: saturation front arrives at well
%% 4 Calculate saturation profile and plots
% 4.1 Calculate travelling distance of saturations
t = linspace(time0, timef, nts);
dfnnt = dfnn(end:-1:index); % inverted sequence of dfnn
Snnt = Snn(end:-1:index); % inverted sequence of Snn

for ti = 1:nts

    for i = 1:(nsn - index + 1)
        Xsf(i, ti) = qt * t(ti) / (A * phi) * dfnnt(i);
    end

end

% 4.2 Plot the derivative curve of fractional flow function
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Derivate of Fractional - Flow Curves:dfw / dSw');
plot(Snn, dfnn, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('derivate of fractional flow function')
xlabel('{\it S}_{nn}');
ylabel('d {\it f}_{nn} / d {\it S}_{nn}');
axis square;
% 4.3 Plot the non-Newtonian fluid saturation profiles
h_fig5 = figure(5);
set(h_fig5, 'color', 'w', 'NumberTitle', 'off', 'Name', 'non - Newtonian Saturation Profiles:Snn(t)');
SwTime = [];

for ti = 1:nts
    plot(Xsf(:, ti), Snnt, '-r', 'LineWidth', 2.0);
    hold on;
    plot([0 Xsf(1, ti)], [1 - Sner 1 - Sner], '-r', 'LineWidth', 2.0);
    SwTime = [SwTime; ['Time = ' num2str(t(1, ti) / 3600, format) ' hour']];
    plot([Xsf(end, ti) Xsf(end, ti)], [Snnf Snnr], '-r', 'LineWidth', 2.0);
    plot([Xsf(end, ti) L], [Snnr Snnr], '-r', 'LineWidth', 2.0);
end

set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('non-Newtonian fluid saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_{nn}');
h_legend5 = legend(SwTime);
set(h_legend5, 'Box', 'on', 'Location', 'best');
axis square;
%% end

%C.2 MATLAB CODE FOR BUCKLEY–LEVERETT SOLUTION:
%BINGHAM NON-NEWTONIAN FLUID

%% Problem description
% Analytical solution of Buckley –Lev erett equation using integral method
% Immiscible displacement of a Bingham non-Newtonian fluid by a Newtonianfluid
% Brooks–Corey relative permeability curve has been applied
close all; clear all; clc;
%% 1 Parameters initialization
% 1.1 rock properties
L = 4.0; % length of domain [m]
A = 1.0; % cross-sectional area [m^2]
phi = 0.20; % porosity of domain
k = 9.869e-13; % permeability of domain [m^2]
alpha = pi * 0.0; % angle between x and horizontal [rad]
% 1.2 fluid properties
rho_ne = 1.0e3; % density of Newtonian fluid [kg/m^3]
rho_nn = 0.9e3; % density of Bingham fluid [kg/m^3]
dlt_rho = rho_ne - rho_nn; % density difference [kg/m^3]
mu_ne = 1.0e-3; % viscosity of Newtonian fluid [Pa*s]
mu_B = 5.0e-3; % Bingham plastic coefficient [Pa*s]
mu_inf = 1.0e20; % infinite viscosity [Pa*s]
G = 5.0e3; % threshold pressure gradient [Pa/m]
% 1.3 relative permeability parameters
Sner = 0.00; % irreducible Newtonian saturation
Snnr = 0.20; % initial non-Newtonian saturation
ne = 2.00; % exponent of Newtonian fluid
nn = 2.00; % exponent of non-Newtonian fluid
krne_max = 0.75; % maximum permeability of Newtonian
krnn_max = 0.75; % maximum permeability of Non-Newtonian
% 1.4 other initial parameters
dlt_Sne = 1e-3; % constant saturation step [-]
Sne = [(Sner + eps):dlt_Sne:(1 - Snnr)]'; % Newtonian fluid saturation vector
g = 9.8067; % gravitational acceleration [m/s^2]
qt = 2e-6; % constant injection rate [m^3/s]
time0 = 86400 * 0.1; % initial calculation time [s]
timef = 86400 * 1.0; % final calculation time [s]
nts = 4; % the number of time steps
nsn = size(Sne, 1); % size of Sne vector
format = '%4.2e'; % precision format for plotting
%% 2 Relative permeabilities, fractional flow function and its derivative
% 2.1 function handles of relative permeability and potential gradient
kr_ne = @(S) krne_max .* ((S - Sner) ./ (1 - Sner - Snnr)) .^ ne; % relative permeability
kr_nn = @(S) krnn_max .* ((1 - S - Snnr) ./ (1 - Sner - Snnr)) .^ nn; % relative permeability
dpx = @(S) -rho_nn * g * sin(alpha) + (qt / A / k + kr_nn(S) / mu_B * G + ... % potential gradient
(kr_ne(S) / mu_ne * rho_ne + kr_nn(S) / mu_B * rho_nn) * g * sin(alpha)) ./ (kr_ne(S) / mu_ne + kr_nn(S) / mu_B);
% 2.2 calculate the non-Newtonian Bingham fluid viscosity
krne = kr_ne(Sne); % relative permeability Newtonian fluid
krnn = kr_nn(Sne); % relative permeability non-Newtonian
px = dpx(Sne); % minus potential gradient
mu_nn = mu_inf * ones(nsn, 1); % initial Bingham fluid viscosity

for i = 1:nsn % Bingham fluid viscosity

    if (abs(px(i)) - G) >= 0.0
        mu_nn(i) = mu_B / (1 - G / abs(px(i)));
    end

end

% 2.3 calculate fractional flow function
mob = (kr_nn(Sne) ./ mu_nn) .* (mu_ne ./ kr_ne(Sne)); % mobility function
fne = 1 ./ (1 + mob) - (A * k .* kr_nn(Sne) ./ (qt .* mu_nn) * ... % fractional flow function
dlt_rho * g * sin(alpha)) ./ (1 + mob);
% 2.4 Evaluate the max Newtonian fluid saturation
for i = 1:nsn

    if fne(i) == 1.0
        Sne_max = Sne(i);
        fne_max = 1.0;
        iSmax = i;
        break;
    end

end

% 2.5 reset the parameters and derivative of fractional flow function
Sne_B = Sne(1:iSmax);
fne_B = fne(1:iSmax);
nsb = size(Sne_B, 1);
dfne_B = ones(nsb, 1);

for i = 3:(nsb - 2)
    dfne_B(i) = (fne_B(i + 1) - fne_B(i - 1)) / (2 * dlt_Sne);
end

dfne_B(2) = (-11 * fne_B(2) + 18 * fne_B(3) - 9 * fne_B(4) + 2 * fne_B(5)) / (6 * dlt_Sne);
dfne_B(1) = abs(2 * dfne_B(2) - dfne_B(3));
dfne_B(nsb - 1) =- (-11 * fne_B(nsb - 1) + 18 * fne_B(nsb - 2) - 9 * fne_B(nsb - 3) + 2 * fne_B(nsb - 4)) / (6 * dlt_Sne);
dfne_B(nsb) = 2 * dfne_B(nsb - 1) - dfne_B(nsb - 2);
%% 3 plot relative permeability, fractional flow and its derivative curves
% 3.1 figure 1: Plot the relative permeability curves
h_fig1 = figure(1);
set(h_fig1, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves: kr');
plot(Sne, krne, '-b', Sne, krnn, '-r', 'LineWidth', 2.0);
hold on;
axis([0.0 1.0 0.0 1.0]);
axis square;
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Brooks-Corey relative permeability curves')
xlabel('{\it S}_{ne}');
ylabel('{\it k}_{r{\it \beta}}');
set(gca, 'YTick', 0:0.2:1);
h_legend1 = legend('Newtonian fluid', 'Bingham fluid');
set(h_legend1, 'Box', 'on', 'Location', 'best');

% 3.2 figure 2: Plot the dp/dx-Sne curves
h_fig2 = figure(2);
set(h_fig2, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Relative Permeability Curves: kr');
plot(Sne, px / 1e5, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('derivates of fractional flow function');
xlabel('{\it S}_{ne}');
ylabel('-d {\it p} / d {\it x}');
axis square;

% 3.3 figure 3: Plot the fractional flow function curve
h_fig3 = figure(3);
set(h_fig3, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Fractional flow function curve');
plot(Sne_B, fne_B, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1], 'YLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('fractional flow function');
xlabel('{\it S}_{ne}');
ylabel('{\it f}_{ne}');
axis square;

% 3.4 Plot the derivate curve of fractional flow function
h_fig4 = figure(4);
set(h_fig4, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Derivatives of Fractional Flow Curves:dfw / dSw');
plot(Sne_B, dfne_B, '-b', 'LineWidth', 2.0);
set(gca, 'XLim', [0 1]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('derivate of fractional flow function')
xlabel('{\it S}_{ne}');
ylabel('d {\it f}_{ne} / d {\it S}_{ne}');
axis square;

%% 4 Calculate advance frontal saturation and travelling distance
dfnet = dfne_B(end:-1:1); % inverted sequence of dfne
fnet = fne_B(end:-1:1); % inverted sequence of fne
Snet = Sne_B(end:-1:1); % inverted sequence of Sne
% 4.1 calculate travelling distance of every saturations
t = linspace(time0, timef, nts);
Sne_t = cell(nts, 1); % saturation Sw,2: initial
Xsf_t = cell(nts, 1); % travelling distance of Sw,2: initial

for ti = 1:nts
    Wi(ti) = qt * t(ti); % injected wetting phase fluid volume

    for i = 1:nsb
        Xsf(i, ti) = Wi(ti) / (A * phi) * dfnet(i); % calculate travelling distance
    end

    Xsf0 = [0; Xsf(1:end - 1, ti)]; % the 2nd travelling distance vector
    dlt_Xsf = Xsf(:, ti) - Xsf0; % dlt_Xsf = Xsf_j - Xsf_(j-1), x_0 = 0
    dlt_Snet = Snet - Sner; % dlt_Sne = Snet,j - Sner

    for j = 1:nsb
        % calculate the injected fluid volume from 0 to Xsf
        Vj = A * phi * sum(dlt_Xsf(1:j) .* dlt_Snet(1:j));

        if Vj >= Wi(ti)
            index(ti) = j; % index of the Snef in Snet
            Snef(ti) = Snet(j); % advance front saturation: Snef
            Xswf(ti) = Xsf(j, ti); % the position of Snef
            break;
        end

    end

    Xsf_t{ti, 1} = Xsf(1:index(ti), ti);
    Sne_t{ti, 1} = Snet(1:index(ti));
end

% 4.2 Evaluate the time of advance frontal saturation Snef
krne_sf = kr_ne(Snef(1));
krnn_sf = kr_nn(Snef(1));
munn_sf = calObjFun(Snef(1), mu_nn, Sne_B);
mob_sf = (krnn_sf / munn_sf) * (mu_ne / krne_sf); % mobility function
fne_sf = 1 ./ (1 + mob_sf) - ... % fractional flow function
A * k .* krnn_sf ./ (qt * munn_sf) * dlt_rho * g * sin(alpha) ./ (1 + mob_sf);
dfne_sf = calObjFun(fne_sf, dfne_B, Sne_B); % derivate of fractional flow function
tf = A * phi * L / (qt * dfne_sf);
%% 5 Newtonian fluid saturation profile and plots
h_fig5 = figure(5);
set(h_fig5, 'color', 'w', 'NumberTitle', 'off', 'Name', 'Newtonian fluid Saturation Profiles:Sne(t)');
SwTime = [];

for ti = 1:nts
    plot(Xsf_t{ti, 1}, Sne_t{ti, 1}, '-r', 'LineWidth', 2.0);
    hold on;
    plot([0 Xsf_t{ti, 1}(1)], [Sne_max Sne_max], '-r', 'LineWidth', 2.0);
    SwTime = [SwTime; ['Time = ' num2str(t(1, ti) / 3600, format) 'hours']];
    plot([Xsf_t{ti, 1}(end) Xsf_t{ti, 1}(end)], [Sne_t{ti, 1}(end) Sner], '-r','LineWidth', 2.0);
    plot([Xsf_t{ti, 1}(end) L], [Sner Sner], '-r', 'LineWidth', 2.0);
end

set(gca, 'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0 L]);
set(gca, 'Fontname', 'Times New Roman', 'FontSize', 10);
title('Newtonian fluid saturation profiles')
xlabel('{\it x} (m)');
ylabel('{\it S}_{ne}');
h_legend2 = legend(SwTime);
set(h_legend2, 'Box', 'on', 'Location', 'best');
axis square;
%% end

%%C.3 AUXILIARY FUNCTIONS

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
            error('Cant find a tangent line through point (Swc, 0). Decrease dlt_Sw !');
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
