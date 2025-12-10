from numpy import array,zeros,abs,max,maximum, minimum









def RungeKutta4_Solver(f, tspan, y0, N):
    
    """
    Parámetros:
        f     : función f(t, y) que define el sistema
        tspan : lista/tupla [t0, tf]
        y0    : condición inicial (array de tamaño 4)
        N     : número de pasos
        
    Retorna:
        t : vector de tiempos
        y : matriz con la solución (N+1, 4)
    """
    
    y = zeros((N+1, 4))             # matriz soluciones
    t = zeros(N+1)                  # vector tiempos
    
    h = (tspan[1] - tspan[0]) / N   # paso
    
    y[0, :] = y0                    # condicion inicial
    t[0] = tspan[0]                 # tiempo inicial
    
    for n in range(N):
        yn = y[n, :]
        tn = t[n]

        k1 = f(tn, yn)
        k2 = f(tn + h/2, yn + h*k1/2)
        k3 = f(tn + h/2, yn + h*k2/2)
        k4 = f(tn + h,   yn + h*k3)

        y[n+1, :] = yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
        t[n+1] = tn + h

    return t, y







    """
    Parámetros:
        f     : función f(t, y) que define el sistema
        tspan : lista/tupla [t0, tf]
        y0    : condición inicial (array de tamaño 4)
        N     : número de pasos
        
    Retorna:
        t : vector de tiempos
        y : matriz con la solución (N+1, 4)
    """

""" def RungeKutta8_Solver(f, tspan, y0, N):


    y = zeros((N+1, 4))
    t = zeros(N+1)

    h = (tspan[1] - tspan[0]) / N

    y[0, :] = y0
    t[0] = tspan[0]

    for n in range(N):
        yn = y[n, :]
        tn = t[n]

        k1  = f(tn + 0,          yn)
        k2  = f(tn + h*(4/27),   yn + (h*(4/27))*k1)
        k3  = f(tn + h*(2/9),    yn + (h/18)*(k1 + 3*k2))
        k4  = f(tn + h*(1/3),    yn + (h/12)*(k1 + 3*k3))
        k5  = f(tn + h*(1/2),    yn + (h/8)*(k1 + 3*k4))
        k6  = f(tn + h*(2/3),    yn + (h/54)*(13*k1 - 27*k3 + 42*k4 + 8*k5))
        k7  = f(tn + h*(1/6),    yn + (h/4320)*(389*k1 - 54*k3 + 966*k4 - 824*k5 + 243*k6))
        k8  = f(tn + h,          yn + (h/20)*(-234*k1 + 81*k3 - 1164*k4 + 656*k5 - 122*k6 + 800*k7))
        k9  = f(tn + h*(5/6),    yn + (h/288)*(-127*k1 + 18*k3 - 678*k4 + 456*k5 - 9*k6 + 576*k7 + 4*k8))
        k10 = f(tn + h,          yn + (h/820)*(1481*k1 - 81*k3 + 7104*k4 - 3376*k5 +
                                              72*k6 - 5040*k7 - 60*k8 + 720*k9))

        y[n+1, :] = yn + (h/840)*(41*k1 + 27*k4 + 272*k5 + 27*k6 + 216*k7 + 216*k9 + 41*k10)
        t[n+1] = tn + h

    return t, y
 """















def RungeKutta45_Solver(f, tspan, y0, N, max_steps = 100000):
    """
    Implementación del método Runge–Kutta–Fehlberg 4(5) (RK45)
    con control de paso adaptativo.

    ----------------------------
    EMBEDDED RUNGE-KUTTA SCHEME RK45
    ----------------------------
    EMBEDDED RUNGE-KUTTA SCHEME RK45 con paso de tiempo adaptable
        k1 = delta_T*F(t_n, dr_n, r_n)
        k2 = delta_T*F(t_n + delta_T/4, dr_n + 1/4*k1, r_n + 1/4*k1)
        k3 = delta_T*F(t_n + delta_T*3/8, dr_n + 3/32*k1 + 9/32*k2, r_n + 3/32*k1 + 9/32*k2)
        k4 = delta_T*F(t_n + delta_T*12/13, dr_n + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3, r_n + 1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3)
        k5 = delta_T*F(t_n + delta_T, dr_n + 439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4, r_n + 439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4)
        k6 = delta_T*F(t_n + delta_T/2, dr_n - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5, r_n - 8/27*k1 + 2*k2 - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5)
        Order 4 estimation--> Uo4_n+1 = U_n + (k1*25/216 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5)
        Order 5 estimation--> Uo5_n+1 = U_n + (k1*16/135 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6)
        Error_array=U_o5_n+1-U_o4_n+1 
        Tol=Tol_a+Tol_r*norm(U_n)
        E = max (Error_array(j)/Tol(j))
        S=(1/E)^(1/5) where 5 is p+1 of the lower order method
        if E<=1: Integration continues
            U_n+1=Uo5_n+1
            delta_T=delta_T*S
        elif 1<=E: Integration is repeated with a lower step
            delta_T=delta_T*S
            i=i-1

    Parámetros:
    f     : función f(t, y) que define el sistema
    tspan : lista/tupla [t0, tf]
    y0    : condición inicial (array de tamaño 4)
    N     : número de pasos   
    
    Retorna:
        t_RK45 : vector de tiempos
        y_RK45 : matriz de soluciones (n_steps, dim)
    """
    

    
    Tol_a = 1e-12                   # Parámetro de tolerancia a
    Tol_r = 1e-12                   # Parámetro de tolerancia r
    
    h = (tspan[1] - tspan[0]) / N   # Paso
    
    y = y0                          # Condicion inicial
    t = tspan[0]                    # Tiempo inicial


    # Guardado de trayectoria
    t_list = [t]
    y_list = [y.copy()]
    
    step_count = 0                  # Se inicializa el contador de pasos

    # Bucle principal
    while t < tspan[1]:                         # Mientras t sea menor al tiempo final

        step_count += 1                         # Incrementa el contador de pasos en 1
        if step_count > max_steps:              # Si se llega al máximo de pasos
            print("RK45: Max steps reached.")   # Mensajito 
            break                               # Se detiene el proceso
        
        # Ajustar último paso
        if t + h > tspan[1]:                    # Si sumarle el paso lo llevaría más lejors del tiempo final         
            h = tspan[1] - t                    # Este se ajusta
        
        E = 2.0                                 # Establece un Error > 1 para entrar al while siguiente siempre


        while E > 1.0:

            # Calcular los k
            k1 = f(t,                   y)
            k2 = f(t + h*1/4,           y + h*(1/4*k1))
            k3 = f(t + h*3/8,           y + h*(3/32*k1 + 9/32*k2))
            k4 = f(t + h*12/13,         y + h*(1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3))
            k5 = f(t + h*1,             y + h*(439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4))
            k6 = f(t + h*1/2,           y + h*(- 8/27*k1 + 2*k2 + - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5))
            
            # Orden 4 (y4) y orden 5 (y5)
            y4 = y + h*(25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5)               # Calcula la solución de orden 4
            y5 = y + h*(16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6) # Calcula la solución de orden 5
            
            # Estimar error
            error = abs(y5 - y4)                            # Array. Calcula el error local tomando la solución de orden 5 como la exacta
            Tol   = Tol_a + Tol_r*maximum(abs(y), abs(y5))  # Array. Maximum devuelve un array [max(1st component),max(2nd component),...]
            E     = max(error / Tol)                        # Scalar. Divide cada componente y selecciona el que tenga el mayor Error_array_j/Tol_j
            
            
            # factor de seguridad
            if E == 0.0:
                S = 2.0                             # Si el error es 0.0, la aproximación es muy buena y S = S_max_value = 2
                                                    # Podría no ponerse esto, porque el máximo será 2 si E es muy pequeño, pero evita probemas de división por 0
            else:
                S = 0.9 * (1.0 / E)**(1.0 / 5.0)    # Si el error no es 0.0, S se calcula
            

            S = maximum(0.1, minimum(S, 2.0))       # El mínimo garantiza que el valor máximo de S sea 2 por si ha subido mucho en el paso anterior. El máximo garantiza que S no sea muy pequeña


            if E > 1.0:

                h = h * S           # Si E > 1, S se utiliza para hacer el tamaño de paso menor
                if h < 1e-15:       # Si el paso se vuelve demasiado pequeño
                    print("El paso temporal en RK45 se ha vuelto demasiado pequeño")
                    t = tspan[1]    # Si h < 1e-15 con E > 1, no converge, se para :(
                    break

            # Si E > 1, el bucle while se repite con un tamaño de paso menor 
        
        t = t + h           # Incrementa el t actual con el delta 
        y = y5              # Guarda la solución de orden 5 como la buena

        y_list.append(y)    # Con .append, guarda todo U para todo t 
        t_list.append(t)    # Con .append, guarda todo t

        h = h * S           # S ajusta el tamaño de paso para el siguiete paso
        




    # ---- Salida ----
    t_RK45 = array(t_list)
    y_RK45 = array(y_list)
    
    return t_RK45, y_RK45





