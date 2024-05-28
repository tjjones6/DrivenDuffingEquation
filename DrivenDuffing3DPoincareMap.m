%% Driven Duffing Oscillator 3D Return Map Animation
% Author:  Tyler Jones
% Contact: tjjones6@wisc.edu
% Date:    05.23.2024

% Description:
% This MATLAB script models the behavior of a driven Duffing oscillator, 
% solving the Driven Duffing Equation using MATLAB's ode45 solver, and
% plotting the resulting trajectory in 3D space with a Poincaré section.

clear all; close all; clc;

%% Parameters
alpha = -1;
beta = 0.25;
gamma = 2.5;
omega = 2;
delta = 0.1;

% Period Definition
T = 2*pi/omega; % Period of the driving force

% Initial Condition (x, x_dot, phi)
initial_condition = [1; 0; 0]; % Initial data: (x=1, x_dot=0, phi=0)

% Time span
time_span = [0 1000]; % Simulation time span

%% Duffing Equation Definition
Duffing_Eq = @(t, Y) [Y(2); 
    gamma*cos(Y(3)) - delta*Y(2) - alpha*Y(1) - beta*Y(1).^3; 
    omega];

% Solve the Duffing Equation using ODE45
[t, Y] = ode45(Duffing_Eq, time_span, initial_condition);

% Axis settings
x_min = min(Y(:,1));
x_max = max(Y(:,1));
y_min = min(Y(:,2));
y_max = max(Y(:,2));
z_min = min(Y(:,3));
z_max = max(Y(:,3));

axis_scalar = 1.25;
x_lim = axis_scalar*[x_min, x_max];
y_lim = axis_scalar*[y_min, y_max];
z_lim = axis_scalar*[z_min, z_max];

%% Plot the 3D Trajectory
figure('units','normalized','Position',[0.1 0.1 .8 .8])
hold on
xlabel('$x$','Interpreter','latex')
ylabel('$\dot{x}$','Interpreter','latex')
zlabel('$\phi$','Interpreter','latex')
title('3D Trajectory of Driven Duffing Oscillator','Interpreter','latex')
grid on; box on; shading interp
camlight HEADLIGHT; lighting gouraud; view([45 45 45])
xlim(x_lim)
ylim(y_lim)
zlim(z_lim)

% Create a plane in the phi = 0 plane
[x_plane, y_plane] = meshgrid(linspace(x_lim(1), x_lim(2), 30), ...
                              linspace(y_lim(1), y_lim(2), 30));
z_plane = zeros(size(x_plane)); % Plane at phi = 0

h_plane = surf(x_plane, y_plane, z_plane, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
colormap turbo

% Plot the trajectory and mark intersections with the Poincaré section
for i = 2:length(Y)
    plot3(Y(1:i,1), Y(1:i,2), Y(1:i,3), 'color', 'b', 'LineWidth', 1.5)
    
    % Check for crossing the phi = 0 plane
    if (Y(i,3) >= 0 && Y(i-1,3) < 0)
        plot3(Y(i,1), Y(i,2), Y(i,3), '*r', 'MarkerSize', 5)
    end
    
    drawnow
end
hold off

%% Plotting Reference Functions
function myfigpref
% MYFIGPREF just makes figures pretty. Written by TGJChandler
% Last edited: 01/01/2018 by TGJChandler
% Comment by Tyler Jones: Thomas Chandler was my professor for math
% 415 (Applied Dynamical Systems, Chaos, and Modeling) @UW-Madison

set(0, 'DefaultAxesFontSize', 20)
set(0, 'DefaultAxesLineWidth', 2);
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultPatchLineWidth', .7);
set(0, 'DefaultLineMarkerSize', 6);

grid on;
box on;

h = gca;
h.TickLabelInterpreter='latex';
h.MinorGridAlpha=0.05;
h.GridAlpha=0.05;
h.FontSize=25;
h.LineWidth=2;

h = gcf;
h.Color = [1,1,1];
end

function fig_xytit(xlab, ylab, tit)
% FIG_XYTIT sets the current figure's xlabel, ylabel, and title in latex format.
% Last edited: 15/06/2021 by TGJChandler

if nargin<3
    tit = '';
end

xlabel(xlab, 'interpreter', 'latex')
ylabel(ylab, 'interpreter', 'latex')
title(tit, 'interpreter', 'latex')
end
