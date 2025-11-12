% Calculo de error, 


clc, clearvars, close all




%% VARIABLES

U0 = [1,0,0,1];     % Vector de estado inicial

T = 200;            % Tiempo total de propagacion en segundos
N = 20000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo


t0 = 0;             % Valor inicial de tiempo


dt1=dt;
dt2= 2*dt;

%% LLAMADA A LAS FUNCIONES

Error(U0, N, dt, t0, 2, @Euler_method);
Error(U0, N, dt, t0, 2, @CrankNicolson_method);
Error(U0, N, dt, t0, 2, @CrankNicolsonMejorado_method, 3);
Error(U0, N, dt, t0, 2, @CrankNicolsonNewton_method, 100, 1e-9);
Error(U0, N, dt, t0, 2, @RungeKutta4_method);
Error(U0, N, dt, t0, 2, @EvilEuler_method);
Error(U0, N, dt, t0, 2, @EvilEulerMejorado_method, 3);


ConvergenceRate(U0, N, dt1, dt2, t0, 2, @Euler_method);
ConvergenceRate(U0, N, dt1, dt2, t0, 2, @CrankNicolson_method);
ConvergenceRate(U0, N, dt1, dt2, t0, 2, @CrankNicolsonMejorado_method, 3);
ConvergenceRate(U0, N, dt1, dt2, t0, 2, @CrankNicolsonNewton_method, 100, 1e-9);
ConvergenceRate(U0, N, dt1, dt2, t0, 2, @RungeKutta4_method);
ConvergenceRate(U0, N, dt1, dt2, t0, 2, @EvilEuler_method);
ConvergenceRate(U0, N, dt1, dt2, t0, 2, @EvilEulerMejorado_method, 3);

%% DECLARACION DE FUNCIONES PRINCIPALES

function Err = Error(U, N, dt, t0, factor, metodo, varargin)
    % metodo: función manejadora (@Euler_method, @CrankNicolson_method, etc.)
    % varargin: argumentos adicionales opcionales (iteraciones, tolerancia, etc.)

    dt2 = dt / factor;   % Paso temporal de la malla fina
    N2 = N * factor;     % Número de pasos de la malla fina

    % Ejecutamos las funciones con los argumentos opcionales (si los hay)
    [U1, pos1, vel1, t1] = metodo(U, N, dt, t0, varargin{:});
    [U2, pos2, vel2, t2] = metodo(U, N2, dt2, t0, varargin{:});

    % Cogemos los valores de la malla fina que coinciden con la gruesa
    U2 = U2(1:factor:end, :);

    q = GiveOrder(metodo);

    % Cálculo del error por Richardson
    Err = dt^q * ((U2 - U1) / (dt^q - dt2^q));

    % Nombre de método legible automáticamente
    nombre = func2str(metodo); 


    % Graficar resultados
    %graficar(pos1, vel1, nombre);
    %graficar_error(Err, t1, nombre);

    % Mostrar el error en consola
    %disp("Error por " + nombre + ":");
    %disp(norm(Err));

end




function q = GiveOrder(method)

    nombre = func2str(method);
    switch nombre
        case 'Euler_method'
            q = 1;
        case {'CrankNicolson_method','CrankNicolsonMejorado_method','CrankNicolsonNewton_method'}
            q = 2;
        case 'RungeKutta4_method'
            q = 4;
        case {'EvilEuler_method','EvilEulerMejorado_method'}
            q = 1;
        otherwise
            q = 1; % por defecto
    end
end


function p = ConvergenceRate(U, N, dt1, dt2, t0, factor, metodo, varargin)

    nombre = func2str(metodo);

    E1 = Error(U, N, dt1, t0, factor, metodo, varargin{:});
    E2 = Error(U, N, dt2, t0, factor, metodo, varargin{:});

    p = log(norm(E1)/norm(E2))/log(dt1/dt2);

    % Mostrar el error en consola
    disp("Ratio de convergencia de " + nombre + ":");
    disp(norm(p));

end


function [U_hist, posiciones, velocidades, t_hist] = Euler_method(U, N, dt, t0)

    t = t0;                     % Le doy al tiempo su valor inicial
    
    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 


    for i = 1:N                             % Inicia el bucle desde el 1 hasta el N
        U = U + dt*F(U, t);                 % Se actualiza el valor de U
        t=t+dt;                             % Se avanza al siguiente instante de tiempo
        t_hist(i+1) = t;                    % Guardamos los valores de tiempo en la variable usada para graficar
        U_hist(i+1, :) = U;                 % Guardamos los valores de la funcion principal en la variable usada para graficar
    end                                     % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

end




function [U_hist, posiciones, velocidades, t_hist] = CrankNicolson_method(U, N, dt, t0)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 
    
    for i = 1:N                                     % Inicia el bucle desde el 1 hasta el N

        U_next_pre = U + dt*F(U, t);                        % Se precalcula el primer valor de Un+1 (se predice con Euler)
        U_next = U + dt/2 * (F(U, t)+F(U_next_pre, t+dt));  % Se hace uso del esquema de CrankNicolson
                                                            % Si se quisiera se podria usar el propio Crank sobre la prediccion varias veces
                                                            % precalculando con el propio Crank y luego ya pasarselo al Crank principal.
                                                            % Esto aumentaria la precision al metodo de CrankNicolson

        t=t+dt;                                     % Se avanza al siguiente instante de tiempo
        U = U_next;                                 % Se iguala U a U_next para el siguiente paso
        t_hist(i+1) = t;                    % Guardamos los valores de tiempo en la variable usada para graficar
        U_hist(i+1, :) = U;                         % Guardamos los valores en la variable usada para graficar

    end                                             % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

end





function [U_hist, posiciones, velocidades, t_hist] = CrankNicolsonMejorado_method(U, N, dt, t0, iter)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           % El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 

    
    for i = 1:N                                                     % Inicia el bucle desde el 1 hasta el N

        U_next_pre = U + dt*F(U, t);                                % Se precalcula el primer valor de Un+1 (se predice con Euler)
        for j = 1:iter
            U_next_pre = U + dt/2 * (F(U, t)+F(U_next_pre, t+dt));  %Se emplea el propio Crank para aproximar mejor el valor de prediccion
        end

        U_next = U + dt/2 * (F(U, t)+F(U_next_pre, t+dt));          % Se hace uso del esquema de CrankNicolson

        t=t+dt;                                     % Se avanza al siguiente instante de tiempo
        U = U_next;                                 % Se iguala U a U_next para el siguiente paso
        t_hist(i+1) = t;                            % Guardamos los valores de tiempo en la variable usada para graficar
        U_hist(i+1, :) = U;                         % Guardamos los valores en la variable usada para graficar

    end                                             % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

end





function [U_hist, posiciones, velocidades, t_hist] = CrankNicolsonNewton_method(U, N, dt, t0, max_iter, tol)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           % El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 

    for i = 1:N
        
        U_pred = U + dt * F(U, t);        % Predicción inicial (Euler) para empezar Newton:

        converged = false;
        for k = 1:max_iter
            G   = U_pred - U - (dt/2) * (F(U, t) + F(U_pred, t + dt));     % Residuo
            JF = JF_Kepler(U_pred);                                         % Jacobiano de F en Z
            JG  = eye(4) - (dt/2) * JF;                                     % Jacobiano de G

            delta = - JG \ G(:);                                            % Paso de Newton: JG * delta = -G

            U_pred_new = U_pred + delta.';                                  % Mantener fila

            if norm(delta) <= tol * (1 + norm(U_pred_new))                  % Chequeo de la tolerancia con absoluto y relativo
                U_pred = U_pred_new;
                converged = true;
                break
            end
            U_pred = U_pred_new;
        end

        if ~converged                   % Si un paso no ha convergido con el maz de iteraciones, se apechuga y se continua
            disp(":(")
        end

        U = U_pred;                     % Acepta U^{n+1}
        t_hist(i+1) = t;                % Guardamos los valores de tiempo en la variable usada para graficar
        U_hist(i+1,:) = U;              % Guarda historial
        t=t+dt;                         % Se avanza al siguiente instante de tiempo
    end

    % Separar para graficar como en tus funciones
    posiciones  = U_hist(:,1:2);
    velocidades = U_hist(:,3:4);
end








function [U_hist, posiciones, velocidades, t_hist] = RungeKutta4_method(U, N, dt, t0)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;     % El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 

    
    for i = 1:N                                     % Inicia el bucle desde el 1 hasta el N

        k1 = F(U, t);                               % Hacemos uso del esquema de RungeKutta 4
        k2 = F(U + dt/2 * k1, t + dt/2);
        k3 = F(U + dt/2 * k2, t + dt/2);
        k4 = F(U + dt * k3, t + dt);

        U = U + dt/6 *(k1 + 2*k2 + 2*k3 + k4);      % Se actualiza el valor de U
        t=t+dt;                                     % Se avanza al siguiente instante de tiempo
        t_hist(i+1) = t;                    % Guardamos los valores de tiempo en la variable usada para graficar
        U_hist(i+1, :) = U;                         % Guardamos los valores en la variable usada para graficar

    end                                             % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

end




function [U_hist, posiciones, velocidades, t_hist] = EvilEuler_method(U, N, dt, t0)


    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           %El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 

    
    for i = 1:N                                 % Inicia el bucle desde el 1 hasta el N

        U_next_pre = U + dt*F(U, t);            % Se precalcula el primer valor de Un+1 (se predice con Euler)
        U_next = U + dt * F(U_next_pre, t+dt);  % Se hace uso del esquema de Euler Inverso
                                                % Si se quisiera se podria usar el propio Crank sobre la prediccion varias veces
                                                % precalculando con el propio Crank y luego ya pasarselo al Crank principal.
                                                % Esto aumentaria la precision al metodo de Euler Inverso

        t=t+dt;                                 % Se avanza al siguiente instante de tiempo
        U = U_next;                             % Se iguala U a U_next para el siguiente paso
        t_hist(i+1) = t;                    % Guardamos los valores de tiempo en la variable usada para graficar
        U_hist(i+1, :) = U;                     % Guardamos los valores en la variable usada para graficar

    end                                         % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

end




function [U_hist, posiciones, velocidades, t_hist] = EvilEulerMejorado_method(U, N, dt, t0, iter)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, 4);     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;     %El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 

    
    for i = 1:N                                     % Inicia el bucle desde el 1 hasta el N

        U_next_pre = U + dt*F(U, t);                    % Se precalcula el primer valor de Un+1 (se predice con Euler)
        for j = 1:iter
            U_next_pre = U + dt * F(U_next_pre, t+dt);  %Se emplea el propio Crank para aproximar mejor el valor de prediccion
        end

        U_next = U + dt * F(U_next_pre, t+dt);          % Se hace uso del esquema de CrankNicolson

        t=t+dt;                                     % Se avanza al siguiente instante de tiempo
        U = U_next;                                 % Se iguala U a U_next para el siguiente paso
        t_hist(i+1) = t;                    % Guardamos los valores de tiempo en la variable usada para graficar
        U_hist(i+1, :) = U;                         % Guardamos los valores en la variable usada para graficar

    end                                             % Termina el bucle


    posiciones = U_hist(:,1:2);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,3:4);    % Del registro guardado cogemos las velocidades

end



























%% DECLARACION DE FUNCIONES AUXILIARES

function F = F(U,t)
    pos = U(1:2);                       % Los dos primeros valores de U son de la posicion
    vel = U(3:4);                       % Los dos ultimos valores de U son de la velocidad

    gamma = 0.3;
    g = 9.81;
    k = 0.5;
    w = 3;

    Ax = 2;
    Ay = 1;

    F = [vel, -gamma*vel(1) - k*pos(1) + Ax*cos(w*t), -g - gamma*vel(2) - k*pos(2) + Ay*sin(w*t)];
    %F = [vel, -2*Ax*wx*vel(1) - wx^2*pos(1), -2*Ay*wy*vel(1) - wy^2*pos(1)];
end

function graficar(posiciones, velocidades, nombre)

    figure;                                                         % Graficar posiciones                                            
    plot(posiciones(:,1), posiciones(:,2), "LineWidth",1.5);        % En el eje x la primera componente, en el eje y la segunda
    xlabel("x");                                                    % Eje x
    ylabel("y");                                                    % Eje y
    title("Trayectoria por "+ nombre);                              % Titulo de la grafica
    grid on;                                                        % Graficamos con un grid puesto
    axis equal;                                                     % Grafica cuadrada

    figure;                                                         % Graficar velocidades                                              
    plot(velocidades(:,1), velocidades(:,2), "LineWidth",1.5);      % En el eje x la primera componente, en el eje y la segunda
    xlabel("x");                                                    % Eje x
    ylabel("y");                                                    % Eje y
    title("Velocidades por "+ nombre);                              % Titulo de la grafica
    grid on;                                                        % Graficamos con un grid puesto
    axis equal;                                                     % Grafica cuadrada

end


function graficar_error(error, t_hist, nombre)


    errNorm = vecnorm(error, 2, 2);   % (N+1 x 1)


    figure;
    %plot(t_hist, errNorm, "r-o","LineWidth",1.5,"MarkerSize",1);
    plot(t_hist, errNorm, "r-o","MarkerSize",1);
    xlabel("Tiempo");
    ylabel("Error");
    title("Estimacion del error de U_h por Richardson en " + nombre);
    grid on;

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