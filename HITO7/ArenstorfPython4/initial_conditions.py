from numpy import array,abs,column_stack
from numpy.linalg import norm


"""
----------------------------
NEWTON SHOOTING METHOD
----------------------------
Problema restringido circular de 3 cuerpos con:
μ = 0.012277471   (Tierra-Luna para la órbita de Arenstorf)
T = 17.0621211718029   (período de la órbita de Arenstorf)

Condiciones iniciales calculadas con el Newton Shooting Method.

- Se realiza una suposición inicial Xk = [x0, y0] (x0 y dy0).
- Se resuelve el problema y, debido a la simetría de la órbita de Arenstorf, en T/2 se debe cumplir que [y, dx] = 0.
- Como la suposición inicial no contiene los valores exactos [x0, dy0], entonces [y, dx] ≠ 0 en T/2.
- El error G en T/2 se evalúa como la diferencia entre la posición exacta y la posición obtenida:  [0, 0] - [y, dx] = -[y, dx].
- El error E es la norma de G.
- Para determinar cuánto debe cambiar Xk, se resuelve:  J * Delta_Xk = G ya que el Jacobiano J mide cómo cambia Xk, y queremos que Xk cambie exactamente G.
- El Jacobiano se evalúa numéricamente mediante diferencias finitas.

"""

def Initial_Conditions(Temporal_Scheme, F, T, N, Initial_guess_x0, Initial_guess_vy0):
    
    Xk_initial_guess = array([Initial_guess_x0,Initial_guess_vy0])    # Guess incial de (x0, dy0/dt)
    Xk = Xk_initial_guess                                             # Guarda el Xk

    Tol = 1e-12
    epsilon = 1e-8
    max_iter = 10000

    U0 = array([Xk[0],0,0,Xk[1]])                                     # Guarda el U0

    T_target = T/2.0                                                  # Integra solo para T/2 por simetría en la que y & dx deben ser = 0
    tspan = array([0,T_target])

    for i in range (max_iter):
        t, U_current = Temporal_Scheme(F, tspan, U0, N)                  # Resuelve para este U0 
        U_current_end = U_current[-1]                                     # Guarda U a T/2
        G = array([U_current_end[1],U_current_end[2]])                    # Guarda U a T/2 en un array G
        Error = norm(G)                                                   # Calcula el error a T/2, ya que las posiciones y & dx deberían ser  =  0 a T/2. entonces, Posiciones = Error

        print(f"Iter {i+1}: x0 = {Xk[0]:.6f}, vy0 = {Xk[1]:.6f} -> Error = {Error:.2e}")

        if Error < Tol: # Si Error < Tol Se ha encontrado U0!!!
            print("Convergencia lograda")
            print(f"U0:")
            print(f"x0 = {U0[0]}")
            print(f"y0 = {U0[1]}")
            print(f"vx0 = {U0[2]}")
            print(f"vy0 = {U0[3]}")
            
            return U0
        
        # Si Error > Tol, U0 no se ha encontrado y el bucle for continua
        # Ahora calculamos el Jacobiano numérico con diferencias finitas
         
        Xk_pert1 = array([Xk[0]+epsilon,Xk[1]])                               # x0 se modifica con epsilon 
        U0_pert1 = array([Xk_pert1[0],0,0,Xk_pert1[1]])                       # nuevo U0 con x0 modificado

        t, U_pert1 = Temporal_Scheme(F,tspan, U0_pert1,N)                    # nuevo U con x0 modificado

        U_pert1_end = U_pert1[-1]                                             # Guarda U a T/2 con x0 modificado
        G_pert1 = array([U_pert1_end[1],U_pert1_end[2]])                      # Calcula el error a T/2 con x0 modificado

        J1 = (G_pert1-G)/epsilon                                              # Calcula la primera columnca del Jacobiano 


        Xk_pert2 = array([Xk[0],Xk[1]+epsilon])                               # dy0 se modifica con epsilon
        U0_pert2 = array([Xk_pert2[0],0,0,Xk_pert2[1]])                       # nuevo U0 con dy0 modificado

        t, U_pert2 = Temporal_Scheme(F,tspan, U0_pert2,N)                    # nuevo U con dy0 modificado

        U_pert2_end = U_pert2[-1]                                             # Guarda U a T/2 con dy0 modificado
        G_pert2 = array([U_pert2_end[1],U_pert2_end[2]])                      # Calcula el error a T/2 con dy0 modificado

        J2 = (G_pert2-G)/epsilon                                              # Calcula la segunda columna del Jacobiano 


        J = column_stack([ J1, J2])         # Construye el Jacobiano numérico
        det_J = J[0,0]*J[1,1]-J[1,0]*J[0,1]   # Calcula el det del Jacobiano

        if abs(det_J) < 1e-15:
            print(f"Iter {i+1}: Error |G| = {Error:.2e}. Singular Jacobiano. Stopping.") # Para si el Jacobiano es singular
            break

        J_inv = array([[J[1,1], -J[0,1]],
                    [-J[1,0], J[0,0]]]) 
        
        J_inv = J_inv/det_J                   # Calcula el inverso del Jacobiano
        delta_X = -J_inv@G                    # Resuleve J_inv*G. delta_Xk mide cuánto debe cambiar Xk
        
        max_step_norm = 0.1                 # Establece un valor máximo para delta_Xk
        current_step_norm = norm(delta_X)

        if current_step_norm > max_step_norm:
            delta_X = delta_X * (max_step_norm / current_step_norm) # Reduce delta_Xk si es mayor que el máximo

        Xk = Xk+delta_X               # Calcula el nuevo x0,dy0
        U0 = array([Xk[0],0,0,Xk[1]]) # Calcula el nuevo U0
        
    return U0








