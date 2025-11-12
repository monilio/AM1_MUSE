

clc, clearvars, close all




%% VARIABLES

U0 = [1, 0];     % Vector de estado inicial

T = 200;            % Tiempo total de propagacion en segundos
N = 2000;           % Numero total de pasos
dt = T/N;           % Paso de tiempo


t0 = 0;             % Valor inicial de tiempo



%{
Euler_method(U0, N, dt, t0);
CrankNicolson_method(U0, N, dt, t0);
CrankNicolsonMejorado_method(U0, N, dt, t0, 3);
CrankNicolsonNewton_method(U0, N, dt, t0, 100, 1e-10);
RungeKutta4_method(U0, N, dt, t0);
EvilEuler_method(U0, N, dt, t0);
EvilEulerMejorado_method(U0, N, dt, t0, 3);
LeapFrog_method(U0, N, dt, t0);
 
%}

region_estabilidad('Euler_method');

region_estabilidad('cranknicolson_method');

region_estabilidad('EvilEuler_method');

region_estabilidad('RungeKutta4_method');


function [U_hist, posiciones, velocidades, t_hist] = Euler_method(U, N, dt, t0)

    t = t0;                     % Le doy al tiempo su valor inicial
    
    U_hist = zeros(N+1, length(U));     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
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


    posiciones = U_hist(:,1);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,2);    % Del registro guardado cogemos las velocidades

    
    graficar(posiciones, velocidades, t_hist, "Euler");
end




function [U_hist, posiciones, velocidades, t_hist] = CrankNicolson_method(U, N, dt, t0)


    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, length(U));     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
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


    posiciones = U_hist(:,1);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,2);    % Del registro guardado cogemos las velocidades

    graficar(posiciones, velocidades, t_hist, 'Crank1');

end





function [U_hist, posiciones, velocidades, t_hist] = CrankNicolsonMejorado_method(U, N, dt, t0, iter)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, length(U));     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
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


    posiciones = U_hist(:,1);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,2);    % Del registro guardado cogemos las velocidades

    graficar(posiciones, velocidades, t_hist, "CrankMejorado");
end





function [U_hist, posiciones, velocidades, t_hist] = CrankNicolsonNewton_method(U, N, dt, t0, max_iter, tol)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, length(U));     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
                                % Esta tiene N+1 filas pues tendra el valor inicial + 1 valor por cada uno de los N pasos
    
    U_hist(1, :) = U;           % El primer valor del historial sera el inicial de U 

    t_hist = zeros(N+1);        % Creo una variable en la que guardar los datos de tiempo cada paso para luego graficarlos
                                
    t_hist(1) = t;              % El primer valor del historial sera el inicial de t 

    for i = 1:N
        
        U_pred = U + dt * F(U, t);        % Predicción inicial (Euler) para empezar Newton:

        converged = false;
        for k = 1:max_iter
            G   = U_pred - U - (dt/2) * (F(U, t) + F(U_pred, t + dt));     % Residuo
            JF = Jacobiano(@(Z) F(Z, t+dt), U_pred);                                         % Jacobiano de F en Z
            JG  = eye(2) - (dt/2) * JF;                                     % Jacobiano de G

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
    posiciones  = U_hist(:,1);
    velocidades = U_hist(:,2);

    graficar(posiciones, velocidades, t_hist, "CrankNewton");
end








function [U_hist, posiciones, velocidades, t_hist] = RungeKutta4_method(U, N, dt, t0)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, length(U));     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
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


    posiciones = U_hist(:,1);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,2);    % Del registro guardado cogemos las velocidades


    graficar(posiciones, velocidades, t_hist, "RK4");
end




function [U_hist, posiciones, velocidades, t_hist] = EvilEuler_method(U, N, dt, t0)


    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, length(U));     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
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


    posiciones = U_hist(:,1);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,2);    % Del registro guardado cogemos las velocidades


    graficar(posiciones, velocidades, t_hist, "EvilEuler");
end




function [U_hist, posiciones, velocidades, t_hist] = EvilEulerMejorado_method(U, N, dt, t0, iter)

    t = t0;                     % Le doy al tiempo su valor inicial

    U_hist = zeros(N+1, length(U));     % Creo una variable en la que guardar los datos de cada paso para luego graficarlos
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


    posiciones = U_hist(:,1);     % Del registro guardado cogemos las posiciones
    velocidades = U_hist(:,2);    % Del registro guardado cogemos las velocidades


    graficar(posiciones, velocidades, t_hist, "EvilEulerMejorado");
end





function [U_hist, posiciones, velocidades, t_hist] = LeapFrog_method(U, N, dt, t0)


    t = t0;                     % Tiempo inicial

    U_hist = zeros(N+1, length(U));     % Historial de estados (posición, velocidad)
    t_hist = zeros(N+1, 1);             % Historial de tiempo

    U_hist(1, :) = U;           % Primer valor: estado inicial
    t_hist(1) = t;              % Primer valor de tiempo

    % Descomponemos U en posición y velocidad inicial
    x = U(1);
    v = U(2);

    % Calcular la aceleración inicial
    a = F(U, t);               % F debe devolver [vx, ax], así que nos interesa la aceleración
    a = a(2);                  % Tomamos la componente de aceleración

    % Paso inicial de medio tiempo para la velocidad
    v_half = v + 0.5 * dt * a;

    % Bucle principal
    for i = 1:N
        % Actualizar posición con la velocidad intermedia
        x = x + dt * v_half;

        % Avanzar el tiempo
        t = t + dt;

        % Calcular la nueva aceleración en el tiempo actualizado
        U_temp = [x, v_half];   % Estado temporal (posicion, velocidad intermedia)
        a_new = F(U_temp, t);
        a_new = a_new(2);

        % Actualizar velocidad intermedia
        v_half = v_half + dt * a_new * 0.5; % medio paso para sincronizar

        % Guardar resultados
        U_hist(i+1, :) = [x, v_half];
        t_hist(i+1) = t;
    end

    % Recuperar posiciones y velocidades para graficar
    posiciones = U_hist(:, 1);
    velocidades = U_hist(:, 2);

    % Graficar resultados
    graficar(posiciones, velocidades, t_hist, "Leap-Frog");
end





















%% DECLARACION DE FUNCIONES AUXILIARES

function F = F(U, t)

    F = [U(2), -U(1)];
end

function graficar(posiciones, velocidades, t_hist, nombre) 
    %--------------------------------------------------------------
    % Graficar posición y velocidad en función del tiempo
    %--------------------------------------------------------------

    figure;                                                        
    plot(t_hist, posiciones, "LineWidth", 1.5);                    
    xlabel("Tiempo t");                                            
    ylabel("Posición x(t)");                                       
    title("Evolución temporal de la posición — " + nombre);        
    grid on;                                                       


%{
     figure;                                                        
    plot(t_hist, velocidades, "LineWidth", 1.5);                   
    xlabel("Tiempo t");                                            
    ylabel("Velocidad v(t)");                                      
    title("Evolución temporal de la velocidad — " + nombre);       
    grid on;                                                       
 
%}

    %--------------------------------------------------------------
    % Diagrama de fase (x vs v)
    %--------------------------------------------------------------

%{
     figure;
    plot(posiciones, velocidades, "LineWidth", 1.5);
    xlabel("Posición x");
    ylabel("Velocidad v");
    title("Diagrama de fase — " + nombre);
    grid on;
    axis equal; 
%}

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





function region_estabilidad(method)
    % REGION_ESTABILIDAD - Grafica la región de estabilidad de un método numérico

    theta = linspace(0, 2*pi, 500);
    r = exp(1i * theta);

    switch lower(method)
        case 'euler_method'
            w = @(r) r - 1;
            name = 'Euler explícito';
        case 'cranknicolson_method'
            w = @(r) (2*(r - 1)) ./ (1 + r);
            name = 'Crank–Nicolson';
        case 'evileuler_method'
            w = @(r) 1 - 1./r;
            name = 'Euler inverso';
        case 'leapfrog_method'
            w = @(r) (r - 1./r) / 2;
            name = 'Leap–Frog';
        case 'rungekutta4_method'
            w = @(r) (r - 1) .* (1 - (r - 1)/2 + (r - 1).^2/6 - (r - 1).^3/24);
            name = 'Runge–Kutta 4';
        otherwise
            error('Método no reconocido: %s', method);
    end

    w_values = w(r);

    figure;
    fill(real(w_values), imag(w_values), 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on;
    plot(real(w_values), imag(w_values), 'k', 'LineWidth', 1.0);
    xlabel('Re(w)');
    ylabel('Im(w)');
    title(['Región de estabilidad: ', name]);
    grid on;
    axis equal;

    % Ejes centrados
    plot(xlim, [0 0], 'k--', 'LineWidth', 1);    % Eje real
    plot([0 0], ylim, 'k--', 'LineWidth', 1);    % Eje imaginario

    % Centrar la vista
    axis([-3 3 -3 3]);   % Ajusta según el método
end