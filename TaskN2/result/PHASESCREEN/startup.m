clear;
close("all");
clc;

% Set up the path for the project
oldpath = path();
oldpath = path(oldpath, ".\src\");
oldpath = path(oldpath, ".\visualize\");
oldpath = path(oldpath, ".\example\");

% Define global constants
c = 299792458.0;
lambda = 1.55e-6;
k = 2*pi/lambda;

% Global plot settings
set(0, "DefaultTextInterpreter", "latex");
set(0, "DefaultAxesFontSize", 12);
% TODO: find a way to set automatically the interpreter of colorbar labels
% set(0, "DefaultColorbarLabelTextInterpreter", "latex");

marker=char('square', 'diamond', 'v', 'x', 'o');
color=char('b', 'r', 'gr', 'm', 'k', 'y');
facecolor=[
    0.749 0.76 0.96;
    0.96 0.77 0.75;
    0.75 0.96 0.78;
    0.88 0.86 0.88;
    0.93 0.75 0.96
];
