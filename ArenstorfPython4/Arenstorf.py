"""
Ampliación de Matemáticas - Master Universitario en Sistemas Espaciales - ETSIAE
Milestone 7: Orbits of the circular restricted three body problem.
1. Integrate Arenstorf’s periodic orbit. Compare results among GBS, RK
   and AM methods.
"""

from numpy import array,concatenate

from initial_conditions import Initial_Conditions
from arenstorf_equations import arenstorf_equations,F_CR3BP

from RK_schemes import RungeKutta45_Solver, RungeKutta4_Solver, RungeKutta8_Solver
from AM_schemes import AdamsMoulton4_Solver
from GBS_schemes import GBS_Solver

from plot_functions import plot_arenstorf
from Richardson_Error import Graficar_Convergencia_Richardson


"""
----------------------------
ARENSTORF ORBIT
----------------------------
Problema restringido de los 3 cuerpos con:
mu = 0.012277471 Tierra - Luna para la Órbita de Arenstorf
T  = 17.0621211718029 Periodo para la Órbita de Arenstorf

"""




#################################################################################
## SELECCIÓN DE CONDICIONES INICIALES
#################################################################################


"""

Las condiciones iniciales de referencia / tabuladas se han obtenido de la siguiente documentación:

https://www.johndcook.com/blog/2020/02/08/arenstorf-orbit/
https://pymgrit.github.io/pymgrit/applications/arenstorf_orbit.html
https://smath.com/ru-RU/file/8YMhEh/GNU-Scientific-Library_-ODE-Solvers_-Arenstorf-orbit_pdf


"""




print("Seleccione la forma de calculo de condiciones iniciales:")
print("1. Obtenidas por NEWTON SHOOTING METHOD")
print("2. TABULADAS")

opcion1 = input("¿Método de obtención de C.I.? (1-2): ")

print("Seleccione la tipología de la orbita:")
print("1. 4 loops")
print("2. 3 loops")
print("3. 2 loops")
opcion2 = input("¿Número de \"loops\" de la órbita? (1-3): ")


match opcion1:
   case "1":
      match opcion2:
         case "1":
            # CI calculada 4 loops
            r0_exact    = array([0.994,0])                              # Definición de la posición inicial de referencia
            v0_exact    = array([0, -2.00158510637908252240537862224])  # Definición de la velocidad inicial de referencia
            U0_exact    = concatenate((r0_exact,v0_exact))              # Condiciones iniciales de referencia
            Initial_guess_x0=0.994
            Initial_guess_vy0=-2.02
            mu          = 0.012277471                                   # mu de Arenstorf
            T           = 17.0652165601579625588917206249               # Periodo de referencia
            N_init      = 1000                                          # Número de pasos para las condiciones iniciales


            U0 = Initial_Conditions(GBS_Solver, arenstorf_equations, T, N_init,Initial_guess_x0,Initial_guess_vy0)

            Error_U0=U0_exact-U0

            print(f"Error a t = 0 respecto a las condiciones iniciales de referencia (x0=0.994, y0=0, vx0=0, vy0=-2.00158510637908252240537862224):")
            print(f"Error de x0 = {Error_U0[0]}")
            print(f"Error de y0 = {Error_U0[1]}")
            print(f"Error de vx0 = {Error_U0[2]}")
            print(f"Error de vy0 = {Error_U0[3]}")
            pass

         case "2":
            # CI calculada 3 loops
            r0_exact    = array([0.994,0])                              # Definición de la posición inicial de referencia
            v0_exact    = array([0, -2.0317326295573368357302057924])   # Definición de la velocidad inicial de referencia
            U0_exact    = concatenate((r0_exact,v0_exact))              # Condiciones iniciales de referencia
            Initial_guess_x0=0.994
            Initial_guess_vy0=-2.0317
            mu          = 0.012277471                                   # mu de Arenstorf
            T           = 11.124340337266085134999734047                # Periodo de referencia
            N_init      = 1000                                          # Número de pasos para las condiciones iniciales


            U0 = Initial_Conditions(GBS_Solver, arenstorf_equations, T, N_init,Initial_guess_x0,Initial_guess_vy0)

            Error_U0=U0_exact-U0

            print(f"Error a t = 0 respecto a las condiciones iniciales de referencia (x0=0.994, y0=0, vx0=0, vy0=-2.00158510637908252240537862224):")
            print(f"Error de x0 = {Error_U0[0]}")
            print(f"Error de y0 = {Error_U0[1]}")
            print(f"Error de vx0 = {Error_U0[2]}")
            print(f"Error de vy0 = {Error_U0[3]}")
            pass
         case "3":
            # CI calculada 2 loops
            r0_exact    = array([1.2,0])                                # Definición de la posición inicial de referencia
            v0_exact    = array([0, -1.049357510])                      # Definición de la velocidad inicial de referencia
            U0_exact    = concatenate((r0_exact,v0_exact))              # Condiciones iniciales de referencia
            Initial_guess_x0=1.2
            Initial_guess_vy0=-1.04
            mu          = 0.012277471                                   # mu de Arenstorf
            T           = 6.192169331                                   # Periodo de referencia
            N_init      = 1000                                          # Número de pasos para las condiciones iniciales


            U0 = Initial_Conditions(GBS_Solver, arenstorf_equations, T, N_init,Initial_guess_x0,Initial_guess_vy0)

            Error_U0=U0_exact-U0

            print(f"Error a t = 0 respecto a las condiciones iniciales de referencia (x0=0.994, y0=0, vx0=0, vy0=-2.00158510637908252240537862224):")
            print(f"Error de x0 = {Error_U0[0]}")
            print(f"Error de y0 = {Error_U0[1]}")
            print(f"Error de vx0 = {Error_U0[2]}")
            print(f"Error de vy0 = {Error_U0[3]}")
            pass

         case _:
            print("Opción no válida.")

   case "2":
      match opcion2:
         case "1":
            # CI de referencia 4 loops
            r0_exact    = array([0.994,0])                              # Definición de la posición inicial de referencia
            v0_exact    = array([0, -2.00158510637908252240537862224])  # Definición de la velocidad inicial de referencia
            U0_exact    = concatenate((r0_exact,v0_exact))              # Condiciones iniciales de referencia
            mu          = 0.012277471                                   # mu de Arenstorf
            T           = 17.0652165601579625588917206249               # Periodo de referencia

            U0 = U0_exact
            pass
         case "2":
            # CI de referencia 3 loops
            r0_exact    = array([0.994,0])                              # Definición de la posición inicial de referencia
            v0_exact    = array([0, -2.0317326295573368357302057924])   # Definición de la velocidad inicial de referencia
            U0_exact    = concatenate((r0_exact,v0_exact))              # Condiciones iniciales de referencia
            mu          = 0.012277471                                   # mu de Arenstorf
            T           = 11.124340337266085134999734047                # Periodo de referencia


            U0 = U0_exact
            pass
         case "3":
            # CI de referencia 2 loops
            r0_exact    = array([1.2,0])                                # Definición de la posición inicial de referencia
            v0_exact    = array([0, -1.049357510])                      # Definición de la velocidad inicial de referencia
            U0_exact    = concatenate((r0_exact,v0_exact))              # Condiciones iniciales de referencia
            mu          = 0.012277471                                   # mu de Arenstorf
            T           = 6.192169331                                   # Periodo de referencia

            U0 = U0_exact
            pass
         case _:
            print("Opción no válida.")
   case _:
      print("Opción no válida.")



N           = 100000                                        # Número de pasos general
num_orbitas = 1                                             # Número de órbitas a realizar
tspan       = array([0,T*num_orbitas]) 





#################################################################################
## PROPAGACIÓN
#################################################################################

t_sol, U_sol = AdamsMoulton4_Solver(arenstorf_equations, tspan, U0, N)
#t_sol, U_sol = GBS_Solver(arenstorf_equations, tspan, U0, N)
#t_sol, U_sol = RungeKutta45_Solver(arenstorf_equations, tspan, U0, N)
#t_sol, U_sol = RungeKutta4_Solver(arenstorf_equations, tspan, U0, N)


U_T=U_sol[-1,0:4]
Error_UT=U0_exact-U_T

print(f"U(T):")
print(f"x(T)={U_T[0]}")
print(f"y(T)={U_T[1]}")
print(f"vx(T)={U_T[2]}")
print(f"vy(T)={U_T[3]}")

print(f"Error a t=T respecto a las condiciones iniciales de referencia (x0=0.994, y0=0, vx0=0, vy0=-2.00158510637908252240537862224):")
print(f"Error de x(T)={Error_UT[0]}")
print(f"Error de y(T)={Error_UT[1]}")
print(f"Error de vx(T)={Error_UT[2]}")
print(f"Error de vy(T)={Error_UT[3]}")






#################################################################################
## CONVERGENCIA POR RICHARDSON
#################################################################################
# Definimos una lista de pasos para ver la evolución.
# Incluimos pasos "malos" (pequeños N) y "buenos" (grandes N)


lista_N = [10,50,100,200,300,500, 1000, 2000, 4000, 5000, 8000,15000,20000,30000,40000]

Graficar_Convergencia_Richardson(
    Temporal_Scheme = AdamsMoulton4_Solver,
   #  Temporal_Scheme = RungeKutta4_Solver,
   F = F_CR3BP,
    U0 = U0,
    Tf_fijo = 1.0,
    lista_N_base = lista_N,
    p_teorico = 4
)

   




#################################################################################
## GRAFICAS
#################################################################################


# ============ 2D =============
plot_arenstorf(U_sol, mu)
