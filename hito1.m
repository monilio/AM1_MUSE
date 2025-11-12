clc, clearvars, close all

% Comienza el primero de los hitos en MATLAB, que servira para demostrar la
% gran superioridad de este lenguje frente a aquellos de los
% revolucionarios. (A menos que me quite demasiado tiempo, en ese caso
% sudare de esta idea, que quiero mucho a mi novia y prefiero aprovechar mi
% tiempo con ella)


%% VARIABLES

U = [1,0,0,1];      % Vector de estado inicial
T = 200;            % Tiempo total de propagacion en segundos
N = 2000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo


U_hist = zeros(N+1, 4);         % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos

U_hist(1, :) = U;               %El primer valor del historial sera el inicial de U 


for i = 1:N                     % Inicia el bucle desde el 1 hasta el N

    pos = U(1:2);                       % Los dos primeros valores de U son de la posicion
    vel = U(3:4);                       % Los dos ultimos valores de U son de la velocidad
    F = [vel, -pos/(norm(pos))^3];      % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    U = U + dt*F;               % Se actualiza el valor de U
    U_hist(i+1, :) = U;         % Guardamos los valores en la variable usada para graficar
end                             % Termina el bucle


posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades


graficar(posiciones, velocidades, "Euler");      % Se llama a la funcion de graficar





U = [1,0,0,1];      % Vector de estado inicial
T = 200;            % Tiempo total de propagacion en segundos
N = 2000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo


U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                            % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos

U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 


for i = 1:N                                     % Inicia el bucle desde el 1 hasta el N

    pos = U(1:2);                                           % Los dos primeros valores de U son de la posicion
    vel = U(3:4);                                           % Los dos ultimos valores de U son de la velocidad
    F = [vel, -pos/(norm(pos))^3];                          % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    U_next_pre = U + dt*F;                                  % Se precalcula el primer valor de Un+1 (se predice con Euler)

    pos_next = U_next_pre(1:2);                             % Los dos primeros valores de U son de la posicion
    vel_next = U_next_pre(3:4);                             % Los dos ultimos valores de U son de la velocidad
    F_next = [vel_next, -pos_next/(norm(pos_next))^3];      % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    U_next = U + dt/2 * (F+F_next);                         % Se hace uso del esquema de CrankNicolson
                                                            % Si se quisiera se podria usar el propio Crank sobre la prediccion varias veces
                                                            % precalculando con el propio Crank y luego ya pasarselo al Crank principal.
                                                            % Esto aumentaria la precision al metodo de CrankNicolson

    U = U_next;                                             % Se iguala U a U_next para el siguiente paso
    U_hist(i+1, :) = U;                                     % Guardamos los valores en la variable usada para graficar

end                                             % Termina el bucle


posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

graficar(posiciones, velocidades, "Crank-Nicolson");      % Se llama a la funcion de graficar






U = [1,0,0,1];      % Vector de estado inicial
T = 200;            % Tiempo total de propagacion en segundos
N = 2000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo

U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                            % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos

U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 



for i = 1:N                                             % Inicia el bucle desde el 1 hasta el N

    pos = U(1:2);                                           % Los dos primeros valores de U son de la posicion
    vel = U(3:4);                                           % Los dos ultimos valores de U son de la velocidad
    F = [vel, -pos/(norm(pos))^3];                          % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    U_next_pre = U + dt*F;                                  % Se precalcula el primer valor de Un+1 (se predice con Euler)

    for j = 1:iter
        pos_next = U_next_pre(1:2);                         % Los dos primeros valores de U son de la posicion
        vel_next = U_next_pre(3:4);                         % Los dos ultimos valores de U son de la velocidad
        F_next = [vel_next, -pos_next/(norm(pos_next))^3];  % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

        U_next_pre = U + dt/2 * (F+F_next);                 % Se hace uso del esquema de CrankNicolson

    end

    pos_next = U_next_pre(1:2);                             % Los dos primeros valores de U son de la posicion
    vel_next = U_next_pre(3:4);                             % Los dos ultimos valores de U son de la velocidad
    F_next = [vel_next, -pos_next/(norm(pos_next))^3];      % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    U_next = U + dt/2 * (F+F_next);                         % Se hace uso del esquema de CrankNicolson
                                                            % Si se quisiera se podria usar el propio Crank sobre la prediccion varias veces
                                                            % precalculando con el propio Crank y luego ya pasarselo al Crank principal.
                                                            % Esto aumentaria la precision al metodo de CrankNicolson



    U = U_next;                                         % Se iguala U a U_next para el siguiente paso
    U_hist(i+1, :) = U;                                 % Guardamos los valores en la variable usada para graficar

end                                                     % Termina el bucle


posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

graficar(posiciones, velocidades, "Crank-Nicolson Mejorado");      % Se llama a la funcion de graficar




U = [1,0,0,1];      % Vector de estado inicial
T = 200;            % Tiempo total de propagacion en segundos
N = 2000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo


U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                            % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos

U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 


for i = 1:N                                     % Inicia el bucle desde el 1 hasta el N

    pos = U(1:2);                               % Los dos primeros valores de U son de la posicion
    vel = U(3:4);                               % Los dos ultimos valores de U son de la velocidad
    F = [vel, -pos/(norm(pos))^3];              % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    k1 = F;                                     % Hacemos uso del esquema de RungeKutta 4


    Uk1 = U + dt/2 * k1;
    pos = Uk1(1:2);                             % Los dos primeros valores de U son de la posicion
    vel = Uk1(3:4);                             % Los dos ultimos valores de U son de la velocidad
    Fk1 = [vel, -pos/(norm(pos))^3];            % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    k2 = Fk1;


    Uk2 = U + dt/2 * k2;
    pos = Uk2(1:2);                             % Los dos primeros valores de U son de la posicion
    vel = Uk2(3:4);                             % Los dos ultimos valores de U son de la velocidad
    Fk2 = [vel, -pos/(norm(pos))^3];            % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    k3 = Fk2;


    Uk3 = U + dt * k3;
    pos = Uk3(1:2);                             % Los dos primeros valores de U son de la posicion
    vel = Uk3(3:4);                             % Los dos ultimos valores de U son de la velocidad
    Fk3 = [vel, -pos/(norm(pos))^3];            % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)

    k4 = Fk3;


    U = U + dt/6 *(k1 + 2*k2 + 2*k3 + k4);      % Actualizamos el valor de U

    U_hist(i+1, :) = U;                         % Guardamos los valores en la variable usada para graficar

end                                             % Termina el bucle


posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

graficar(posiciones, velocidades, "Runge-Kutta 4");      % Se llama a la funcion de graficar




%% DECLARACION DE FUNCIONES AUXILIARES

function graficar(posiciones, velocidades, nombre)

    figure;                                                         % Graficar posiciones                                            
    plot(posiciones(:,1), posiciones(:,2), 'LineWidth',1.5);        % En el eje x la primera componente, en el eje y la segunda
    xlabel('x');                                                    % Eje x
    ylabel('y');                                                    % Eje y
    title('Trayectoria por '+ nombre);                              % Titulo de la grafica
    grid on;                                                        % Graficamos con un grid puesto
    axis equal;                                                     % Grafica cuadrada

    figure;                                                         % Graficar velocidades                                              
    plot(velocidades(:,1), velocidades(:,2), 'LineWidth',1.5);      % En el eje x la primera componente, en el eje y la segunda
    xlabel('x');                                                    % Eje x
    ylabel('y');                                                    % Eje y
    title('Velocidades por '+ nombre);                              % Titulo de la grafica
    grid on;                                                        % Graficamos con un grid puesto
    axis equal;                                                     % Grafica cuadrada

end


