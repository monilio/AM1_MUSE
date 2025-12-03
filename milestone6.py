#Usados en el trabajo final de AM1


from numpy import array,concatenate,zeros,abs,max,log,real,imag,isclose,eye,block,all,meshgrid,linspace,finfo,copy,sqrt, maximum, minimum
from numpy.linalg import norm,eigvals
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import inspect
from differential_equations import N_body_problem
from temporal_schemes import Euler, RK4, Crank_Nicolson, Evil_Euler, Leapfrog
from aux_functions import Cauchy_problem
from numpy import linspace, array



def RungeKutta45_Solver(f, y0, tspan, N=100, max_steps = 100000):


    
    Tol_a = 1e-12   # Parámetro de tolerancia a
    Tol_r = 1e-12   # Parámetro de tolerancia r
    
    
    h = (tspan[1] - tspan[0]) / N   # paso
    
    y = y0                    # condicion inicial
    t = tspan[0]                 # tiempo inicial


    # Guardado de trayectoria
    t_list = [t]
    y_list = [y.copy()]
    
    step_count = 0
    # Bucle principal
    while t < tspan[1]:

        step_count += 1 # Increases step count by 1
        if step_count > max_steps:
            print("RK45: Max steps reached.")
            break
        
        # Ajustar último paso
        if t + h > tspan[1]:
            h = tspan[1] - t
        

        E = 2.0 # Establishes an Error > 1 to enter the next while in every step


        while E > 1.0:

            # Calcular los k
            k1 = f(t,                   y)
            k2 = f(t + h*1/4,           y + h*(1/4*k1))
            k3 = f(t + h*3/8,           y + h*(3/32*k1 + 9/32*k2))
            k4 = f(t + h*12/13,         y + h*(1932/2197*k1 - 7200/2197*k2 + 7296/2197*k3))
            k5 = f(t + h*1,             y + h*(439/216*k1 - 8*k2 + 3680/513*k3 - 845/4104*k4))
            k6 = f(t + h*1/2,           y + h*(- 8/27*k1 + 2*k2 + - 3544/2565*k3 + 1859/4104*k4 - 11/40*k5))
            
            # Orden 4 (y4) y orden 5 (y5)
            y4 = y + h*(25/216*k1 + 1408/2565*k3 + 2197/4104*k4 - 1/5*k5)               # Calculates order 4 solution
            y5 = y + h*(16/135*k1 + 6656/12825*k3 + 28561/56430*k4 - 9/50*k5 + 2/55*k6) # Calculates order 5 solution
            
            # Estimar error
            error = abs(y5 - y4)                            # Array. Calculates the local error considering the order 5 solution as the exact one
            Tol   = Tol_a + Tol_r*maximum(abs(y), abs(y5))  # Array. Maximum returns an array [max(1st component),max(2nd component),...]
            E     = max(error / Tol)                        # Scalar. Divides each component and selects the one in which Error_array_j/Tol_j is bigger
            
            
            # factor de seguridad
            if E == 0.0:
                S = 2.0                             # If error is 0.0, the aproximation is very good and S = S_max_value = 2
            else:
                S = 0.9 * (1.0 / E)**(1.0 / 5.0)    # If error is not 0.0, S is calculated
            

            S = maximum(0.1, minimum(S, 2.0))       # The minimum guarantees S_max_value is 2 as it could be bigger in the previous calculation. Maximum guarantees S is not too small


            if E > 1.0:

                h = h * S           # If E > 1, S is used to make the time step smaller
                if h < 1e-15:
                    print("El paso temporal en RK45 se ha vuelto demasiado pequeño")
                    t = tspan[1]    # If h < 1e-15 with E > 1, no convergence, stop
                    break

            # If E > 1, the while repeats with a smaller step 
        
        t = t + h           # Increases actual t by delta_T
        y = y5              # Saves the order 5 solution as actual U

        y_list.append(y)    # With .append, it saves ALL U for all t 
        t_list.append(t)    # With .append, it saves ALL t

        h = h * S           # S adjusts the step size for the new step
        




    # ---- Salida ----
    t_RK45 = array(t_list)
    y_RK45 = array(y_list)
    
    return t_RK45, y_RK45






def CR3BP_equations(t, vec, mu = 0.012277471):
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
d2x/dt=dx/dt=0
d2x/dt=dx/dt=0

Then: x -(1-mu)*(x+mu)/sqrt((x+mu)^2+y^2)^3 -mu*(x-1+mu)/sqrt((x-1+mu)^2+y^2)^3=0
      y -(1-mu)*y/sqrt((x+mu)^2+y^2)^3 -mu*y/sqrt((x-1+mu)^2+y^2)^3=0
"""
def Lagrange_Points(mu):
    
    def Equation(x):
        return x - (1 - mu) * (x + mu) / (x + mu)**3 - mu * (x - 1 + mu) / (x - 1 + mu)**3
    x_L1_guess = 1 - mu - 0.1 
    L1_x = fsolve(Equation, x_L1_guess)[0]
    x_L2_guess = 1 - mu + 0.1 
    L2_x = fsolve(Equation, x_L2_guess)[0]
    x_L3_guess = -mu - 1.05 
    L3_x = fsolve(Equation, x_L3_guess)[0]
    
    x_L45 = 0.5 - mu
    y_L4 = sqrt(3) / 2
    y_L5 = -sqrt(3) / 2
    
    L_points = {
        'L1': (L1_x, 0.0),
        'L2': (L2_x, 0.0),
        'L3': (L3_x, 0.0),
        'L4': (x_L45, y_L4),
        'L5': (x_L45, y_L5),
    }
    
    return L_points




r0 = array([0.487,0.86]) 
dr0 = array([0, 0]) 
delta_T=0.001 
N=10000
U0=concatenate((r0,dr0)) 
mu=0.01215058
T = 15
tspan = array([0,T])


L_points=Lagrange_Points(mu)
L1_x, L1_y = L_points['L1']
L2_x, L2_y = L_points['L2']
L3_x, L3_y = L_points['L3']
L4_x, L4_y = L_points['L4']
L5_x, L5_y = L_points['L5']

print(f"L1: x = {L1_x:.6f}, y = 0.000000")
print(f"L2: x = {L2_x:.6f}, y = 0.000000")
print(f"L3: x = {L3_x:.6f}, y = 0.000000")
print(f"L4: x = {L4_x:.6f}, y = {L4_y:.6f}")
print(f"L5: x = {L5_x:.6f}, y = {L5_y:.6f}")



t_sol, U_sol = RungeKutta45_Solver(CR3BP_equations, tspan, U0, N)
