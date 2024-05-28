%% Driven Duffing Oscillator (Varying Parameter)
% Author:  Tyler Jones
% Contact: tjjones6@wisc.edu
% Date:    04.09.2024 

%{
This MATLAB script models the behavior of a driven Duffing oscillator, a 
nonlinear dynamical system common in engineering and physics. By numerically 
solving the Driven Duffing Equation using MATLAB's ode45 solver, it 
generates a Poincaré section plot, illustrating the oscillator's phase space 
dynamics. 

See the following for more info: 
<https://en.wikipedia.org/wiki/Duffing_equation> 
Steven H. Strogatz: Nonlinear Dynamics and Chaos
%}

clear all; close all; clc;

%% Parameter Library and User Input
%{ 
Driven Duffing Oscillator Equation:
$\ddot{x} + \delta\dot{x} + \alpha x + \beta x^{3}=\gamma \cos(\omega t)$

alpha: Coefficient of linear damping/stiffness
beta:  Coefficient of cubic damping/stiffness
gamma: Amplitude of external driving force
omega: Angular frequency of external driving force
delta: Coefficient of velocity damping/stiffness
%}

% Initial Condition (x,x_dot,phi)
y0 = [1; 0]; % Initial data/condition: (x=1, x_dot=0, phi=0)

% Axis Controls
b = 7;
axis_controls = [-b b];
bound = b;

%% Vector Field Generation
%{
This section computes the vector field for the driven Duffing oscillator, 
defining functions for the first-order differential equations of position 
and velocity. It then generates a grid of x and y values and evaluates the 
velocity components.
%}

%% Simulation and Visualization Loop
%{
This section initializes a figure for visualization and sets up parameters 
for the simulation, including defining the Duffing equation using anonymous
functions and configuring video writing settings. It iterates through 
increasing values of the phase angle `phi`, and plotting the resulting 
Poincaré section.

System Definition:
\dot{x} = v
\dot{v} = \gamma\cos(\omega) - \delta\dot{x} - \alpha x - \beta x^3
\dot{\phi} = \omega
%}

figure('units','normalized','Position',[0.1 0.1 .8 .8])

% myWriter = VideoWriter('DuffingEquation3.mp4', 'MPEG-4');
% myWriter.FrameRate = 60;
% open(myWriter);

max_time = 500;
gamma_final = 2.5;
c = 0.1;
tspan = [0,max_time];
for i = 1:ceil(gamma_final/c)
    clf

    alpha = -1;
    beta = 1;
    gamma = c*i;
    omega = 1.4;
    delta = 0.1;

    % Solve Duffing Equation via ODE45
    Duff_Eq = @(t, Y) [Y(2); gamma*cos(omega*t) - delta*Y(2) - alpha*Y(1) - beta*Y(1).^3];
    
    % ODE45 Solver
    [t, Y] = ode45(Duff_Eq, tspan, y0);
    
    % Plot phase space trajectory
    plot(Y(:,1), Y(:,2),'LineWidth',1); % x vs x_dot
    hold on
    plot(Y(1,1), Y(1,2), 'k*'); % Initial point
    grid on
    fig_xytit('$x$','$\dot{x}$')
    title(['Phase Space: $\ddot{x} + \delta\dot{x} + \alpha x + \beta x^3 = \gamma \cos(\omega t)$, $\gamma = $ ',num2str(gamma)],'Interpreter','latex');
    xlim(axis_controls)
    ylim(axis_controls)
    hold off
    pause(0.01)

%     frame = getframe(gcf);
%     writeVideo(myWriter,frame);
end
% close(myWriter)

%% Plotting Reference Functions
function myfigpref
%   MYFIGPREF just makes figures pretty. Written by TGJChandler
%
%   Last edited: 01/01/2018 by TGJChandler
%
%   Comment by Tyler Jones: Thomas Chandler was my professor for math
%   415 (Applied Dynamical Systems, Chaos, and Modeling) @UW-Madison

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
% FIG_XYTIT sets the current figure's xlabel, ylabel, and title
%   in latex format.
%
%   Last edited: 15/06/2021 by TGJChandler

if nargin<3
    tit = '';
end

xlabel(xlab,'interpreter','latex')
ylabel(ylab,'interpreter','latex')
title(tit,'interpreter','latex')
end
