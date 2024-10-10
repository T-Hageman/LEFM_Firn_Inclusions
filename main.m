clearvars
close all
clc

if (isempty(gcp('nocreate')))
    parpool('Threads')
end

%% Material paramters for user to change
E1 = 9.5;       % Young's modulus of consolidated glacial ice E_i [GPa]
E0 = 8;         % Difference in Young's modulus of consolidated glacial ice and Young's modulus of firn E_i-E_s [GPa]
D = 32.5;       % Tuned constant for property distribution [m]
H = 125;        % Glacier thickness [m]
v = 0.35;       % Poisson ratio [-]
rho_i = 917;    % Density of fully consolidated ice [kg/m^3]
rho_sea = 1020; % Density of oceanwater [kg/m^3]
hw = 0.5.*H;    % Oceanwater height [m]
rho_s = 350;    % Density of firn layer [kg/m^3]
g = 9.81;       % Gravitational acceleration [m/s^2]
kcrit = 1.00e5; % Critical stress intensity factor []
x = 0.0;        % Meltwater depth ratio [-]
rho_w = 1000;   % Density of meltwater [kg/m^3]

%% Perform LEFM code for varying meltwater depth ratios and glacier heights
resolution = 20; % resolution of the plots (number of line segments)

for H = 125 % Heights for which to plot results
    parfor x1 = 1:resolution+1 % index for plotting
        x = (x1-1)/(resolution); % Parameter sweep of meltwater depth ratio
        
        L = LEFM_Density(v, D, H, rho_i, rho_sea, hw, rho_s, g, kcrit, x);
        L1 = LEFM_Isotropic(v, H, rho_i, rho_sea, hw, g, kcrit, x);
        L2 = LEFM_Modulus(E1, E0, D, v, H, rho_i, rho_sea, hw, g, kcrit, x);
        L3 = LEFM_Density_Modulus(E1, E0, D, v, H, rho_i, rho_sea, rho_s, hw, g, kcrit, x);
        
        X = [' x = ', num2str(x), ' H = ', num2str(H) , ' d/H = ',  num2str(L), ' d1/H = ',  num2str(L1), ' d2/H = ',  num2str(L2), ' d3/H = ',  num2str(L3)];
        disp(X)
        
        MeltwaterDepth(x1) = x;
        CrevDepth_Density(x1,1)         = L;  % Crevasse depth data outputted for plotting - VARIABLE DENSITY
        CrevDepth_Isotropic(x1,1)       = L1; % Crevasse depth data outputted for plotting - ISOTROPIC ICE
        CrevDepth_Modulus(x1,1)         = L2; % Crevasse depth data outputted for plotting - VARIABLE MODULUS
        CrevDepth_DensityModulus(x1,1)  = L3; % Crevasse depth data outputted for plotting - VARIABLE DENSITY & MODULUS
    end
    PlotResults(MeltwaterDepth, CrevDepth_Density, CrevDepth_Isotropic, CrevDepth_Modulus, CrevDepth_DensityModulus)
    title('H='+string(H)+" m", 'FontSize',18,'Interpreter','latex');
end

%% Figure plotting
function PlotResults(MeltwaterDepth, CrevDepth_Density, CrevDepth_Isotropic, CrevDepth_Modulus, CrevDepth_DensityModulus)
    figure()
    
    box on
    set(gca,'FontSize',12)
    set(gca, 'FontName', 'Times New Roman')
    set(gca, 'LabelFontSizeMultiplier', 1.5)
    set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
    
    plot(MeltwaterDepth, CrevDepth_Density, '-r','LineWidth',1, 'DisplayName', 'H = 125 m $\rho(z)$');
    hold on,  
    plot(MeltwaterDepth, CrevDepth_Isotropic, '-b','LineWidth',1, 'DisplayName', 'H = 125 m $\rho$ = 917 kg m$^{-3}$ E = 9.5 GPa');
    plot(MeltwaterDepth, CrevDepth_Modulus, '-k','LineWidth',1, 'DisplayName', 'H = 125 m $E(z)$');
    plot(MeltwaterDepth, CrevDepth_DensityModulus, '-g','LineWidth',1, 'DisplayName', 'H = 125 m $\rho(z)$ $E(z)$');
    
    xlim([-0.05 1.05])
    ylim([-0.05 1.05])
    xticks(0:0.2:1)
    yticks(0:0.1:1)
    xlabel('Meltwater Depth Ratio ($h_\mathrm{s}$/$d_\mathrm{s}$)', 'FontSize',18,'interpreter','latex')
    ylabel('Normalised Crevasse Depth($d_\mathrm{s}$/$H$)', 'FontSize',18,'interpreter','latex')
    
    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.XLabel.FontSize = 18;
    ax.YAxis.FontSize = 10;
    ax.YLabel.FontSize = 18;

    legend('FontSize',10,'interpreter','latex','Location','southeast');
end
