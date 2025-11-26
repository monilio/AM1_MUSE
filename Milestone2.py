from numpy import array, concatenate, zeros, linspace
from numpy.linalg import norm, solve, LinAlgError
import matplotlib.pyplot as plt
from aux_functions import  Newton



def Euler(F, U, t1, t2, **kwargs):

    return U + (t2-t1)*F(U, t1)





def Crank_Nicolson(F, U, t1, t2, tol_J=None, max_iter=None, tol_N=None, **kwargs):

    def crank_nicolson_function(U_next, F, U_n, t1, t2):                # funcion a resolver del Crank-Nicolson

        dt = t2-t1
        return (0.5*(F(U_next, t2) + F(U_n, t1))) - ((U_next - U_n)/dt) # 0.5*(F(U^(n+1),t(n+1)) + F(U^(n),t(n))) - (U^(n+1) - U^(n))/dt = 0

    sol = Newton(crank_nicolson_function, U, tol_J, max_iter, tol_N, F, U, t1, t2) # Se resuelve con un Newton y palante

    return sol





def RK4(F, U, t1, t2, **kwargs):    # El esquema de toda la vida
    dt = t2 - t1
    k1 = F(U, t1)
    k2 = F(U+0.5*k1*dt, t1+0.5*dt)
    k3 = F(U+0.5*k2*dt, t1+0.5*dt)
    k4 = F(U+k3*dt, t1+dt)

    return U + (1.0/6.0)*dt*(k1 + 2*k2+2*k3 + k4)






def Evil_Euler(F, U, t1, t2, tol_J=None, max_iter=None, tol_N=None, **kwargs):


    def Evil_euler_function(U_next, F, U, t1, t2):  # Funcion que devuelve la funcion a resolver del Euler inverso
        dt = t2 - t1
        return U_next-U-dt*F(U_next, t2)            # U^(n+1) - U^(n) - dt*F(U^(n+1), t(n+1)) = 0

    sol = Newton(Evil_euler_function, U, tol_J, max_iter, tol_N, F, U, t1, t2) # Se resuelve con un Newton y palante

    return sol






def Cauchy_problem(F, U0, t, temporal_scheme, **kwargs):

    N = len(t)
    N_v = len(U0)
    U = zeros((N, N_v))
    U[0,:] = U0

    for n in range(0, N-1):
        U[n+1,:] = temporal_scheme(F, U[n,:], t[n], t[n+1], **kwargs)

    return U





def F(U, t):

    r = U[0:2]
    rd = U[2:4]

    return concatenate((rd, -r/norm(r)**3), axis=None)





def plotter(U, method, N, T, dt):
    plt.axis("equal")
    plt.plot(U[:, 0], U[:, 1], label=f"N={N}, dt={dt}")
    plt.scatter(U[:, 0], U[:, 1], color='teal', s=2, label="Time steps")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(method)
    plt.show()
    return

U0 = array([1, 0, 0, 1])
T = 10
N = 20000
t = linspace(0, T, N)
dt = (t[1]-t[0])/N
U = Cauchy_problem(F, U0, t, Evil_Euler, tol_J=1e-10, max_iter=11111, tol_N=1e-10)

plotter(U, "Evil Euler", N, T, dt)