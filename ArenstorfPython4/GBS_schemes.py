from numpy import array, inf

def GBS_Solver(f, tspan, y, N, tol=1e-11):


    """
    ============================================================================
    INFO
        Esta función como tal no calcula nada, solo ahce el bucle que permite ir actualizando las cosas
        Las que de verdad calculan cosas son GBS_Step y GBS_Midpoint
        Estas se llaman dentro de los buscles de esta función

    INPUTS
        f       es la función F de Arenstorf
        tspan   es [t0, tf] (tiempos inicial y final)
        y       es la condición inicial (vector fila)
        N       es el número de pasos inicial 
        tol     es la tolerancia deseada (opcional, por defecto 1e-11)
    ============================================================================
    """

    
    h = (tspan[1]- tspan[0])/N  # Paso inicial. Este cambiará porque el GBS tiene un paso de tiempo dinámico.
    y_GBS = [y]                 # Primer valor de la solución (la hacemos de tamaño variable, por eso este es solo el primer valor)
    t_GBS = [tspan[0]]          # Primer valor de tiempos (t0) (lo hacemos de tamaño variable, por eso este es solo el primer valor)

    m_list = array([2, 4, 6, 8, 10, 12, 14, 16])    # Lista de subdivisiones (para saber más, mirar GBS_Step y GBS_Midpoint)

    max_iter = 100000000   # Límite para evitar bucles infinitos
    t = tspan[0]           # Instante de tiempo del paso de integración, comienza con el valor inicial


    while t < tspan[1] and len(t_GBS) < max_iter:   # Mientras el tiempo no supere al tiempo final
                                                    # Y mientras la longitud del vector de tiempos no supere al valor de tieraciones máximas

        # El paso es dinámico, pero hay que hacer que se ajuste de tal manera que el último paso acabe justo en el tiempo final
        if t + h > tspan[1]:    # Si el tiempo actual más el paso que se va a dar suepra al tiempo final       
            h = tspan[1] - t    # Se cambia el tamaño del paso para que acabe justo en el tiempo final


        success, y_new, h_new = GBS_Step(f, t, y, h, m_list, tol)       # Se llama a la función que calcula la solución y el paso usado para hallarla
                                                                        # Si el paso sí era lo suficientemente bueno, devuelve succes = true
                                                                        # Si el paso no era lo suficientemente bueno, devuelve succes = false
                                                                        # Si el paso ha sido rechazado, también da un h_new diferente al introducido
        if success:             # Si el paso ha sido aceptado

            t = t + h           # Se actualiza el valor del tiempo
            y = y_new           # Se actualiza el valor de vector de estado
            t_GBS.append(t)     # Se introduce el valor de tiempo en nuestro vector dinámico de tiempos
            y_GBS.append(y)     # Se introduce el valor del vector de estado en nuestro vector dinámico de soluciones
            h = h_new           # Se actualiza el valor del propio tamaño del paso para la siguiente iteración

        else:                   # Si el paso ha sido rechazado

            h = h_new           # Se actualiza también el valor del tamaño de paso pero no se da el paso como tal
                                # Esto hace que en la siguiente iteración se pruebe este nuevo paso para ver si es suficiente o no
                                # Y así ocurre un bucle hasta que el h_new sea lo suficientemente bueno

    y_GBS = array(y_GBS)
    return t_GBS, y_GBS









def GBS_Step(f, tn, yn, h, m_list, tol=1e-11):

    """
    ============================================================================
    INFO
        Esta función da el paso avanzando al próximo valor del vector estado
        Aun así, hace uso de otra función, llamada GBS_Midpoint

    INPUTS
        f       es la función F de Arenstorf
        tn      es el tiempo actual
        yn      es el vector estado actual
        h       es el tamaño de paso inicial
        m_list  es la lista de subdivisiones
        tol     es la tolerancia deseada 
    ============================================================================
    """

    
    kmax = len(m_list)          # Se toma la longitud de m_list como el tamaño de la matriz de extrapolación triangular
    
    Extrap = [[None for _ in range(kmax)] for _ in range(kmax)]     # Se genera la tabla de extrapolación triangular con el tamaño especificado
                                                                    # Esta está inicialmente llena únicamente de ceros
                                                                    # Extrap(i,j): fila i, columna j
                                                                    # Esta será una matriz en la que cada elemento será un vector de estado

    success = False       # Por ahora no se ha logrado el paso, hasta que se diga lo contrario
    y_best = yn.copy()    # Se guarda en y_best el valor del vector de estado de partida en este paso

    
    for i in range(kmax):                       # Se hace una iteración por cada subdivisión 

        m = m_list[i]                           # Se toma el valor correspondiente en la lista de subdivisiones

        Y_mid = GBS_Midpoint(f, tn, yn, h, m)   # Se calcula el vector estado de punto medio necesitado llamando a una función
                                                # Esto hace el paso integrando con m subpasos desde t hasta t + h
                                                # Claramente usando un subpaso de tamaño h/m
                                                
        Extrap[i][0] = Y_mid                    # Este vector estado se guarda en la priemera de las columnas de la fila en la que se esté en la mtriz de extrapolación
                                                # Esta primera columna, siempre tendrá un vector de estado simple para el que no se ha usado extrapolación



        # Aquí se realiza la extrapolación de Richardson
        for j in range(1, i+1):

            """
            ============================================================================
            Básicamente es la extrapolación de Richardson de toda la vida. Es una forma de calcular una solución "y"
            haciendo uso de y(t+h) e y(t+h/m). De esta manera, para cada valor de "m" se obtendrán un valores diferentes 
            en la extrapolación. Es por eso que para cada "m" hay una fila. Es decir, "m" es como el orden de extrapolación
            porque determina el tamaño de la subdivisión de h

            Lo chulo de la matriz es que en vez de hacer la extrapolación entre y(t+h) e y (t+h/m) cambiando solo el "m",
            haces la extrapolación entre y(t+h/(m_i)) e y(t+h(m_i+1)). Es decir la extrapolación siempre se hace entre 
            ese valor de "m" y el de la fila anterior a esta para dar la columna siguiente.

            De esta manera, se tiene que cada fila tiene un número de columnas no nulo igual a "m" (matriz triangular).

            El orden del error de de cada extrapolación es del orden de 2*i

            Fórmula de extrapolación:
            # E(i,j) = E(i,j-1) + (E(i,j-1) - E(i-1,j-1)) / ( (m_list(i)/m_list(i-j+1))^2 - 1 )
            ============================================================================
            """

            factor = (m_list[i] / m_list[i - j])**2 - 1
            Extrap[i][j] = Extrap[i][j-1] + (Extrap[i][j-1] - Extrap[i-1][j-1]) / factor

        



        # Aquí se realiza la estimación del error
        if i > 0:                                           # Si el el orden de "m" en el que se está no es 1
            err = max(abs(Extrap[i][i] - Extrap[i][i-1]))   # Se calcula el error entre la extrapolación de la última columna y la penúltima
        else:                                               # Si estás en el primer orden de "m" y solo tienes una columna (m = 1)
            err = inf                                       # Se pone que el error es infinito
        



        # Aquí se hace una comprobación de tolerancia
        if err < tol:               # Si este error calculado es menor que la tolerancia dada                                  
            success = True          # Entonces el paso ha sido un éxito
            y_best = Extrap[i][i]   # Se toma que el valor de y es el de la última extrapolación de esta fila
            p = 2 * (i + 1)                 
            #h_new = h * min(1.001, max(0.0001, (tol/err)^(1/(p+1))))
            if err == 0:
                h_new = h * 1.001      # Evita división por cero
            else:
                h_new = h * min(1.001, (tol / err)**(1/(p+1)))

            """
            ============================================================================
            Esta última línea es importante.
            Esta  actualiza el valor del tamaño del paso de tiempo en función del error obtenido.

            la función de (tol/err)^(1/(p+1)) básicamente hace que vayas a un tamaño de paso h que haga que el próximo error que nos
            de, sea lo más similar posible a la tolerancia. Esto hace que siempre vayas al límite, tratando de conseguir el mayor 
            tamaño de paso que cumpla la tolerancia. 

            Esto aumenta el paso siempre, pero lo hace de menor medida cuanto más cerca hayas estado de la tolerancia.
            Para no liarla y que el paso aumente una barbaridad, se tiene cuidado y se pone un máximo. Si no, si un paso 
            te fuera muy bien, podría dar un siguiente paso enorme que desastibilizase todo. 
            ============================================================================
            """

            return success, y_best, h_new
      

        # NOTA: Nótese que el bucle no siempre llega al último valor de "m", si consigue estar por debajo de la tolerancia
        #       con un "m" menor, no es necesario continuar con "m"s más grandes
  

    h_new = h * 0.25    # Si incluso con el mayor "m" no se ha logrado estar por debajo de la tolerancia, se actualiza el tamaño del paso de tiempo
                        # Este nuevo valor se pasa a la función principal para que desde este momento, siempre se parta desde este
                        # Esto posiblemente sea optimizable haciendo uso de (tol/err)^(1/(p+1))

    
    return success, y_best, h_new











def GBS_Midpoint(f, tn, yn, h, m):

    """
    ============================================================================
    INFO
        Esta función calcula el Midpoint, que será el valor de la primera columna
        de la matriz de extrapoación

    INPUTS
        f       es la función F de Arenstorf
        tn      es el tiempo actual
        yn      es el vector estado actual
        h       es el tamaño de paso inicial
        m       es el valor de subdivisión de tamaño de paso h
    ============================================================================
    """

    delta = h/m                     # Tamaño del subpaso
    Y_minus1 = yn.copy()            # Valor de "y" inicial

    Y0 = yn + delta * f(tn, yn)     # Valor de lo que sería el primer paso de "y" usando un tamaño de paso h/m con Euler
                                    
    """
        k1 = f(tn, yn)
        k2 = f(tn + delta/2, yn + delta*k1/2)
        k3 = f(tn + delta/2, yn + delta*k2/2)
        k4 = f(tn + delta,   yn + delta*k3)

        Y0 = (yn + (delta/6)*(k1 + 2*k2 + 2*k3 + k4)) 
    """

                                

    # Aquí se hace un leap-frog
    for k in range(2, m+1):                    # Una iteración por cada una de las "m" subdivisiones
        t_k = tn + (k-1)*delta                 # Tiempo t_{k}
        Y1 = Y_minus1 + 2*delta * f(t_k, Y0)   # Se actualiza la solución dando un doble paso de tiempo con la función evaluada en el punto medio temporal de los tiempos (t) y (t+2*h/m))
        Y_minus1 = Y0                          # Se actualiza para que este sea el paso anterior (de k-1)
        Y0 = Y1                                # Se actualiza para que sea el nuevo paso medio

        """
        ============================================================================
        Para entender mejor esto, sobretodo esas dos ultimas actualizaciones, se hace un esquema: 
        O ------ | ------- O 

        El primer O sería Y_minus1, es decir tu paso anterior (k-2)
        El | sería Y0, ese punto medio Y(k-1) que usas entre (t) y (t+2*delta), es decir, la función evaluada tras un delta (h/m) respecto del tiempo de Y(k-2)
        El segundo O establece el Y1, solución nueva a la que llegas Y(k)

        Así en el sisguinte paso, con Y_minus1 = Y0 y Y0 = Y1 acabas dando este salto:
        O ------- | ------- O 
                  O ------- | ------- O 

        Y se repite
        El primer O sería Y_minus1, es decir tu paso anterior (k-2)
        El | sería Y0, ese punto medio Y(k-1) que usas entre (t) y (t+2*delta), es decir, la función evaluada tras un delta (h/m) respecto del tiempo de Y(k-2)
        El segundo O establece el Y1, solución nueva a la que llegas Y(k)
        Así en el sisguinte paso, con Y_minus1 = Y0 y Y0 = Y1 acabas dando este salto:
        O ------- | ------- O 
                  O ------- | ------- O 
                            O ------- | ------- O 
        
        ...
        ============================================================================
        """
    

    t_final = tn + h   # Tras esto se actualiza el tiempo con el paso completo

    Y_end = 0.5 * (Y0 + Y_minus1 + delta * f(t_final, Y0))  # Hace referencia al "modified midpoint"
                                                            # Básicamente hace una corrección de Gragg al final

    return Y_end

