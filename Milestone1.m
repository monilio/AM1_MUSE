clc, clearvars, close all

% Comienza el primero de los hitos en MATLAB, que servira para demostrar la
% gran superioridad de este lenguaje frente a aquellos de los
% revolucionarios. (A menos que me quite demasiado tiempo, porque tengo que cotizar)


%% VARIABLES

U = [1,0,0,1];      % Vector de estado inicial
T = 200;            % Tiempo total de propagacion en segundos
N = 2000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo



%% LLAMADA A LAS FUNCIONES


Euler_method(U, N, dt);
CrankNicolson_method(U, N, dt);
CrankNicolsonMejorado_method(U, N, dt, 3);
RungeKutta4_method(U, N, dt);
CrankNicolsonNewton_method(U, N, dt, 100, 1e-10)



%% DECLARACION DE FUNCIONES PRINCIPALES

function Euler_method(U, N, dt)


    U_hist = zeros(N+1, 4);         % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                    % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;               %El primer valor del historial sera el inicial de U 


    for i = 1:N                     % Inicia el bucle desde el 1 hasta el N
        U = U + dt*F(U);            % Se actualiza el valor de U
        U_hist(i+1, :) = U;         % Guardamos los valores en la variable usada para graficar
    end                             % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades


    graficar(posiciones, velocidades, "Euler");      % Se llama a la funcion de graficar
end





function CrankNicolson_method(U, N, dt)

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 

    
    for i = 1:N                                     % Inicia el bucle desde el 1 hasta el N

        U_next_pre = U + dt*F(U);                   % Se precalcula el primer valor de Un+1 (se predice con Euler)
        U_next = U + dt/2 * (F(U)+F(U_next_pre));   % Se hace uso del esquema de CrankNicolson
                                                    % Si se quisiera se podria usar el propio Crank sobre la prediccion varias veces
                                                    % precalculando con el propio Crank y luego ya pasarselo al Crank principal.
                                                    % Esto aumentaria la precision al metodo de CrankNicolson

        U = U_next;                                 % Se iguala U a U_next para el siguiente paso
        U_hist(i+1, :) = U;                         % Guardamos los valores en la variable usada para graficar

    end                                             % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

    graficar(posiciones, velocidades, "Crank-Nicolson");      % Se llama a la funcion de graficar
end




function CrankNicolsonMejorado_method(U, N, dt, iter)

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 

    
    for i = 1:N                                             % Inicia el bucle desde el 1 hasta el N

        U_next_pre = U + dt*F(U);                           % Se precalcula el primer valor de Un+1 (se predice con Euler)
        for j = 1:iter
            U_next_pre = U + dt/2 * (F(U)+F(U_next_pre));   %Se emplea el propio Crank para aproximar mejor el valor de prediccion
        end
        U_next = U + dt/2 * (F(U)+F(U_next_pre));           % Se hace uso del esquema de CrankNicolson


        U = U_next;                                         % Se iguala U a U_next para el siguiente paso
        U_hist(i+1, :) = U;                                 % Guardamos los valores en la variable usada para graficar

    end                                                     % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

    graficar(posiciones, velocidades, "Crank-Nicolson Mejorado");      % Se llama a la funcion de graficar
end




function CrankNicolsonNewton_method(U, N, dt, max_iter, tol)


    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 

    for i = 1:N
        
        U_pred = U + dt * F(U);        % Predicción inicial (Euler) para empezar Newton:

        converged = false;
        for k = 1:max_iter
            G   = U_pred - U - (dt/2) * (F(U) + F(U_pred));     % Residuo
            %JF = JF_Kepler(U_pred);                             % Jacobiano de F en Z
            JF = Jacobiano(@(Z) F(Z), U_pred);
            JG  = eye(4) - (dt/2) * JF;                         % Jacobiano de G

            delta = - JG \ G(:);                                % Paso de Newton: JG * delta = -G

            U_pred_new = U_pred + delta.';                      % Mantener fila

            if norm(delta) <= tol * (1 + norm(U_pred_new))      % Chequeo de la tolerancia con absoluto y relativo
                U_pred = U_pred_new;
                converged = true;
                break
            end
            U_pred = U_pred_new;
        end

        if ~converged       % Si un paso no ha convergido con el maz de iteraciones, se apechuga y se continua
            disp(":(")
        end

        U = U_pred;                    % Acepta U^{n+1}
        U_hist(i+1,:) = U;             % Guarda historial
    end

    % Separar para graficar como en tus funciones
    posiciones  = U_hist(:,1:2);
    velocidades = U_hist(:,3:4);
    graficar(posiciones, velocidades, "Crank–Nicolson (Newton)");
end




function RungeKutta4_method(U, N, dt)

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 

    
    for i = 1:N                                     % Inicia el bucle desde el 1 hasta el N

        k1 = F(U);                                  % Hacemos uso del esquema de RungeKutta 4
        k2 = F(U + dt/2 * k1);
        k3 = F(U + dt/2 * k2);
        k4 = F(U + dt * k3);

        U = U + dt/6 *(k1 + 2*k2 + 2*k3 + k4);      % Actualizamos el valor de U

        U_hist(i+1, :) = U;                         % Guardamos los valores en la variable usada para graficar

    end                                             % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

    graficar(posiciones, velocidades, "Runge-Kutta 4");      % Se llama a la funcion de graficar
end



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


function F = F(U)
    pos = U(1:2);                       % Los dos primeros valores de U son de la posicion
    vel = U(3:4);                       % Los dos ultimos valores de U son de la velocidad
    F = [vel, -pos/(norm(pos))^3];      % F = (x_dot, -x/(pos1^2 + pos2^2)^1.5)
end


function J = JF_Kepler(U)   % Jacobiano de F(U) para Kepler 2D:
                            % F(U) = [vx; vy; ax; ay],  a = -r / ||r||^3  con r = [x;y]

    x = U(1);
    y = U(2);
    r2 = x^2 + y^2;
    r  = sqrt(r2);
    if r < 1e-12
        error('Singularidad: ||pos|| ~ 0');     % Yatusabe por si las moscas
    end
    r3 = r^3;
    r5 = r^5;

    a_xx = -1/r3 + 3*x*x/r5;        %df3(U)/dx
    a_xy = 3*x*y/r5;
    a_yx = a_xy;
    a_yy = -1/r3 + 3*y*y/r5;

    J = [ 0   0   1   0
          0   0   0   1
          a_xx a_xy 0  0
          a_yx a_yy 0  0 ];
end


function J = Jacobiano(F, U)

    % Evalúa F en el punto base
    F0 = F(U);
    n = length(U);
    m = length(F0);
    J = zeros(m, n);

    h = 1e-8; % Paso de perturbación
    for j = 1:n
        Upert = U;

        % Pequeña perturbación relativa
        if abs(U(j)) > 1e-8
            h_j = h * abs(U(j));
        else
            h_j = h;
        end

        Upert(j) = Upert(j) + h_j;
        F1 = F(Upert);
        J(:, j) = (F1 - F0) / h_j;
    end
end