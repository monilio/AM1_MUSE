from numpy import array,concatenate,zeros,sqrt,log, linalg, array,log10
import matplotlib.pyplot as plt


#####################################################
##### LAS DOS USAN LAS MISMAS ECUACIONES
#####################################################

def arenstorf_equations(t, vec, mu = 0.012277471):
    """
    Ecuaciones del problema restringido de 3 cuerpos (CR3BP)
    en el marco rotatorio (órbita de Arenstorf).
    
    Parámetros:
        t   : tiempo (no se usa, pero se mantiene para compatibilidad)
        vec : vector de estado [x, y, x_dot, y_dot]
        mu  : parámetro de masa

    Retorna:
        dvec_dt : derivadas del estado
    """

    x, y, x_dot, y_dot = vec

    D1 = ((x + mu)**2 + y**2)**1.5
    D2 = ((x - 1 + mu)**2 + y**2)**1.5

    x_dot_dot = x + 2*y_dot - (1 - mu)*(x + mu)/D1 - mu*(x - 1 + mu)/D2
    y_dot_dot = y - 2*x_dot - (1 - mu)*y/D1      - mu*y/D2

    sol = array([x_dot, y_dot, x_dot_dot, y_dot_dot])

    return sol





"""
----------------------------
PROBLEMA RESTRINGIDO CIRCULAR DE 3 CUERPOS
----------------------------
Problema general de N cuerpos:
d^2r_i/dt^2 = - sum(G * m_j * (r_i - r_j) / |r_i - r_j|^3)
  (la suma es sobre j distinto de i, hasta N)

Consideramos dos masas grandes M1 y M2, y una tercera muy pequeña m3 ≈ 0.  
Con esta suposición, M1 y M2 describen una órbita kepleriana de dos cuerpos.

Unidades adimensionales:
 - M1 + M2 = 1
 - La distancia entre M1 y M2 es R = 1
 - La constante gravitatoria G = 1
 - La velocidad angular (para la órbita circular) w = 1

Parámetro de masa: mu = M2 / (M1 + M2) → M1 = 1 - mu, M2 = mu

El sistema de referencia rota con M1 y M2 con velocidad angular w.  
M1 y M2 permanecen fijos en este sistema de referencia:  
M1 = (-mu, 0),  M2 = (1 - mu, 0).  
El origen coincide con el centro de masas.

El problema estudiará el movimiento de la masa m3. La ecuación es:

d^2r/dt^2 = -∇(pot) - 2 w * dr/dt  
donde -2 w * dr/dt es la fuerza de Coriolis.  
En el caso circular:
   Fx = 2 * dy/dt  
   Fy = -2 * dx/dt

y el potencial es:  
pot = pot_grav + pot_centrífugo  
    = -GM1 / sqrt((x - x1)^2 + y^2) -GM2 / sqrt((x - x2)^2 + y^2) - 1/2 (x^2 + y^2)

Sustituyendo posiciones:  
pot = -(1 - mu) / sqrt((x + mu)^2 + y^2) - mu / sqrt((x - 1 + mu)^2 + y^2) - 1/2 (x^2 + y^2)

Entonces:

d^2x/dt^2 =  2 * dy/dt - d(pot)/dx  
        = 2 * dy/dt + x - (1 - mu)(x + mu) / [( (x + mu)^2 + y^2 )^(3/2)] - mu (x - 1 + mu) / [( (x - 1 + mu)^2 + y^2 )^(3/2)]

d^2y/dt^2 = -2 * dx/dt - d(pot)/dy  
        = -2 * dx/dt + y - (1 - mu) y / [( (x + mu)^2 + y^2 )^(3/2)] - mu y / [( (x - 1 + mu)^2 + y^2 )^(3/2)]
"""
def F_CR3BP(t, Y):
    """
    Adaptación del Problema de los 3 Cuerpos Restringido Circular.
    Entrada:
        t: Tiempo (aunque es autónomo, los solvers lo requieren)
        Y: Vector de estado [x, y, vx, vy]
    Salida:
        dY: Derivada del estado [vx, vy, ax, ay]
    """
    # Definir mu (puedes cambiarlo según el sistema, ej. Tierra-Luna)
    mu = 0.012277471 
    
    # Desempaquetar Y en posición (r) y velocidad (dr)
    r = Y[0:2]   # r[0] es x, r[1] es y
    dr = Y[2:4]  # dr[0] es vx, dr[1] es vy

    # Distancias a los primarios (para no repetir el cálculo del sqrt**3)
    d1 = ((r[0] + mu)**2 + r[1]**2)**(1.5)
    d2 = ((r[0] - 1 + mu)**2 + r[1]**2)**(1.5)

    # F1 es simplemente la velocidad (dr)
    F1 = dr
    
    # F2 es la aceleración (tus ecuaciones)
    F2 = zeros(2)
    F2[0] = 2*dr[1] + r[0] - (1-mu)*(r[0]+mu)/d1 - mu*(r[0]-1+mu)/d2
    F2[1] = -2*dr[0] + r[1] - (1-mu)*r[1]/d1 - mu*r[1]/d2
    
    # Concatenar para devolver un vector de tamaño 4
    F_resultante = concatenate((F1, F2))
    
    return F_resultante











def  df_dY(vec, mu = 0.012277471):

    """
    Derivada manual de la orbita de arenstorf
    
    Parámetros:
        vec : vector de estado [x, y, x_dot, y_dot]

    Retorna:
        df_dY : derivada manual de las ecuaciones
    """

    # Variables
    x = vec[0]
    y = vec[1]
    xdot = vec[2]
    ydot = vec[3]

    # Distancias
    r1sq = (x + mu)**2 + y**2
    r2sq = (x - 1 + mu)**2 + y**2

    r1 = sqrt(r1sq)
    r2 = sqrt(r2sq)

    D1 = r1**3
    D2 = r2**3

    # Factores para derivadas
    r1_5 = r1**5
    r2_5 = r2**5

    # ------------------------------
    # Inicializar Jacobiano
    # ------------------------------
    Jf = zeros((4, 4))

    # ========================================
    # 1) x' = xdot
    # ========================================
    Jf[0,0] = 0
    Jf[0,1] = 0
    Jf[0,2] = 1
    Jf[0,3] = 0

    # ========================================
    # 2) y' = ydot
    # ========================================
    Jf[1,0] = 0
    Jf[1,1] = 0
    Jf[1,2] = 0
    Jf[1,3] = 1

    # ========================================
    # 3) xdot'
    # ========================================

    # Parciales de la parte gravitatoria para r1
    Ax1 = ( -2*(x+mu)**2 + y**2 ) / r1_5      # d/dx ((x+mu)/r1**3)
    Ay1 = ( -3*y*(x+mu) ) / r1_5            # d/dy ((x+mu)/r1**3)

    # Parciales de la parte gravitatoria para r2
    Ax2 = ( -2*(x-1+mu)**2 + y**2 ) / r2_5    # d/dx term r2
    Ay2 = ( -3*y*(x-1+mu) ) / r2_5

    # dx_dot' / d(x,y,xdot,ydot)
    Jf[2,0] = 1 - (1-mu)*Ax1 - mu*Ax2
    Jf[2,1] =     - (1-mu)*Ay1 - mu*Ay2
    Jf[2,2] = 0
    Jf[2,3] = 2

    # ========================================
    # 4) ydot'
    # ========================================

    # d/dx (y/r1**3)
    Bx1 = -3*y*(x+mu)/r1_5

    # d/dy (y/r1**3)
    By1 = (r1sq - 3*y**2) / r1_5

    # d/dx (y/r2**3)
    Bx2 = -3*y*(x-1+mu)/r2_5

    # d/dy (y/r2**3)
    By2 = (r2sq - 3*y**2) / r2_5

    # dy_dot' / d(x,y,xdot,ydot)
    Jf[3,0] =     - (1-mu)*Bx1 - mu*Bx2
    Jf[3,1] = 1 - (1-mu)*By1 - mu*By2
    Jf[3,2] = -2
    Jf[3,3] = 0

    return Jf