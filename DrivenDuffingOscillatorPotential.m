%% Driven Duffing Oscillator (Potential/Energy Surface)
% Author:  Tyler Jones
% Contact: tjjones6@wisc.edu
% Date:    04.010.2024

clear all; close all; clc

% Parameter Library
alpha = -1;
beta = 0.25;
gamma = 2.5;
omega = 2;
delta = 0.1;

%% Simulation and Visualization Loop

t_final = 10;
x_dot = @(x,v) v;
v_dot = @(x,v,t) gamma*cos(omega*t) - delta*v - alpha*x - beta*x.^3;

% ---
% Generate mesh
% ---
x = linspace(-10,10,200);
v = linspace(-10,10,200);
t = linspace(0,t_final,100);
[X,V] = meshgrid(x,v);
UU = x_dot(X,V);

% Begin time stepping
figure('units','normalized','Position',[0.1 0.1 .8 .8])

myWriter = VideoWriter('DuffingEquationPotential.mp4', 'MPEG-4');
myWriter.FrameRate = 30;
open(myWriter);

for i = 1:length(t)
    clf;

    energy = gamma*cos(omega*i*t_final/length(t)) - delta*V - alpha*X - beta*X.^3;
    VV = v_dot(X,V,i*t_final/length(t));

    % Plot the potential surface
    subplot(1,2,1)
    hold on
    myfigpref
    fig_xytit('$x$','$v$','Potential Surface')
    meshc(X,V,energy)
    view([131 60 45])
    colorbar
    colormap jet
    contour(X,V,energy,[0 0],'LineColor','white','LineWidth',5);
    axis auto
    hold off

    subplot(1,2,2)
    hold on
    myfigpref
    fig_xytit('$x$','$v$','Potential Contour')
    contourf(X,V,energy)
    strmslice1 = streamslice(X,V,UU,VV,5);
    set(strmslice1, 'Color', 'white','linewidth',0.5);
    colorbar
    colormap jet
    shading interp
    contour(X,V,energy,[0 0],'LineColor','white','LineWidth',5);
    axis auto
    hold off

    pause(0.01)

    frame = getframe(gcf);
    writeVideo(myWriter,frame);

end

close(myWriter)

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
