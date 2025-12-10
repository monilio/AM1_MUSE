from numpy import zeros, transpose, inf, linalg, identity
from numpy.linalg import norm
from arenstorf_equations import df_dY

def AdamsMoulton4_Solver(f, tspan, y0, N):

    """
    ============================================================================
    INFO
        Esta función presenta el algoritmo de Admas Moulton 4. Este necesita determinar unos pasos iniciales y luego
        llamar a Newton para resolver una ecuación implícita. La ecuación que busca resolver este método es:

        y(n+1) = y(n) + (h/24) * [9*f(t(n+1), y(n+1)) + 19*f(t(n), y(n)) - 5*f(t(n-1), y(n-1)) + f(t(n-2), y(n-2))]

        Como se puede ver, se necesitan los valores de f en n, n-1 y n-2 para avanzar a n+1

        Es de orden 4 y requiere 3 valores previos

    INPUTS
        f       es la función F de Arenstorf
        tspan   es [t0, tf] (tiempos inicial y final)
        y       es la condición inicial (vector fila)
        N       es el número de pasos inicial 
    ============================================================================
    """

    y_AM4 = zeros((N+1, 4))         # Vector de soluciones
    t_AM4 = zeros(N+1)              # Vector de tiempos

    h = (tspan[1] - tspan[0]) / N   # Paso de timepo

    y_AM4[0] = y0                   # Condición inicial
    t_AM4[0] = tspan[0]             # Tiempo inicial



    # Los pasos iniciales que se necesitan se pueden calcular con Euler
    """
     for n = 1:2                                    
        yn = y_AM4(n,:) 
        tn = t_AM4(n) 

        y_AM4(n+1,:)      =   yn + h * f(tn, yn) 
        t_AM4(n+1)    =   tn + h 
    end 
    """

    for n in range(2):                                                 #  Mejor que un Euler, se puede usar un RK4
        yn = y_AM4[n]
        tn = t_AM4[n]

        k1 = f(tn, yn)
        k2 = f(tn + h/2, yn + h*k1/2)
        k3 = f(tn + h/2, yn + h*k2/2)
        k4 = f(tn + h,   yn + h*k3) 

        y_AM4[n+1] = yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4)    # Se calcula hasta n = 3, para así tener n-1 y n-2
        t_AM4[n+1] = tn + h                    
    


    # Aquí se hace el bucle principal
    for n in range(2, N): 

        # Se establecen los valores previos n-1 y n-2
        yn = y_AM4[n, :] 
        tn = t_AM4[n]

        tn_p1 = tn + h 
        tn_m1 = tn - h 
        tn_m2 = tn - 2*h 

        y_m1 = y_AM4[n-1,:] 
        y_m2 = y_AM4[n-2,:]     


        # Se plantea la ecuación de Adams-Moulton que deberá resolver Newton
        # F = @(Y) Y - yn - (h/24) * (9 * f(tn_p1, Y) + 19 * f(tn,yn) - 5 * f(tn_m1,y_m1) + 1 * f(tn_m2,y_m2)) 
        def F(Y):
            return Y - yn - (h/24)*(
                9  * f(tn_p1, Y) +
                19 * f(tn,    yn) -
                5  * f(tn_m1, y_m1) +
                     f(tn_m2, y_m2)
            )


        # Se le da un paso inicial a Newton calculado mediante un Euler explícito
        #Y0 = yn + h * f(tn, yn) 

        k1 = f(tn, yn)                                      # Meh, vamos a probar con un RK4 también
        k2 = f(tn + h/2, yn + h*k1/2) 
        k3 = f(tn + h/2, yn + h*k2/2) 
        k4 = f(tn + h,   yn + h*k3) 
        Y0 = (yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4)) 
        

        
        # Se llama a Newton para que lo resuelva
        yn_p1 = Newton_Solver(F, Y0, h, J_flag = True)  # SI le paso "true" usa el jacobiano manual


        # Guardar solución
        y_AM4[n+1,:] = yn_p1 
        t_AM4[n+1] = tn_p1 

    return t_AM4, y_AM4














def Newton_Solver(F, Y0, h, tol = 1e-12, maxIter = 2000, J_flag = False):

    """
    ============================================================================
    INFO
        Esta función resuelve ecuaciones implícitas mediante el método de Newton.

        Newton busca resolver la ecuación F(Y) = 0 de manera iterativa, usando esta ecuación: 

        Y_{k+1} = Y_k - J(Y_k)^{-1} * F(Y_k)

        Donde J es el Jacobiano de F.
        

    INPUTS
        F       es la función F de Arenstorf
        Y0      es el tiempo actual
        tol     es la tolerancia deseada 
        maxIter es el número máximo de iteraciones para no entrar en bucles demasiado largos
    ============================================================================
    """



    Y = Y0              # Ponerlo así hace que sea un vector columna

    n = len(Y)          # Número de variables del sistema (dimensión del vector Y)

    hJ = 1e-11          # Paso para diferenciar

    for k in range(maxIter):        # Se realizan una serie de ireaciones hasta que se cumpla el máximo o se cumpla la tolerancia
        Fval = F(transpose(Y))      # Se evalúa la función F en el punto Y actual
        
        if norm(Fval, inf) < tol:   # Si el módulo máximo de Fval es menor que la tolerancia    
            return Y                # Se ha convergido, se sale de la función  
                                    # Se recuerda que lo que se busca resolver el Fval = 0, por eso este es el criterio de convergencia 
                                 

        # Aquí se realiza el jacobiano numérico mediante la definición de la derivada

        if J_flag:                                          # Si se ha indicado que se use el jacobiano manual
            J = identity(n) - (3*h/8) * df_dY(transpose(Y)) # Se usa el jacobiano manual calculado con la función df_dY()
                                                            # Esta es la forma del jacobiano para las ecuaciones de Arenstorf

        else:                                   # Si no se ha indicado que se use el jacobiano manual
            J = zeros((n, n))                   # Se inicializa la matriz jacobiana       
            for j in range(n):                  # Bucle para cada columna del jacobiano
                Yp = Y.copy()
                Yp[j] = Yp[j] + hJ
                J[:,j] = transpose(F(transpose(Yp)) - Fval) / hJ # Derivada numérica por definición
            
        

        try:
            dY = linalg.solve(J, -transpose(Fval))
        except linalg.LinAlgError:
            raise RuntimeError("El Jacobiano es singular en Newton.")
        

        
        Y = Y + dY         # Actualizar el valor de Y

        
        if norm(dY, inf) < tol: # Esta también vale como condición de parada
                                # Si no gustase, se pueden poner tolerancias diferentes para las dos condiciones de parada
            return Y            # Se ha convergido, se sale de la función
        
    

    print('Newton no convergió en las iteraciones máximas especificadas.')   # Si no se logra converger en el número máximo de iteraciones, salta un aviso
    return Y