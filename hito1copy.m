clc, clearvars, close all



%% VARIABLES

U = [1,0,0,1];      % Vector de estado inicial
T = 200;            % Tiempo total de propagacion en segundos
N = 2000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo



pos = U(1:2); % x y
vel = U(3:4); % x_dot y_dot

F = [vel, -pos/norm(pos)^3];

for i = 1:N




end

