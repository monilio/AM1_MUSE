from numpy import sqrt, real, imag, array, linspace
from numpy.random import uniform
from numpy import meshgrid
from temporal_schemes import Euler, Crank_Nicolson, RK4, Evil_Euler, Leapfrog
from aux_functions import Cauchy_problem
from differential_equations import oscilador
import matplotlib.pyplot as plt

Puntos_graficar = [array([0.5,0]), array([1,0]),  array([2,0]), array([3,0]), array([4, 0])]
T = 10
N = 1000
t = linspace(0, T, N)
tol_jacobian = 1e-9
N_max = 10000
newton_tol = 1e-10


def plot_orbit(problem, Puntos_graficar, t, time_scheme, **kwargs):

    for U0 in Puntos_graficar:
        U = Cauchy_problem(problem, U0, t, time_scheme, **kwargs)
        plt.plot(U[:,0], U[:,1], marker="o", color='red')
        plt.text(x=2.3, y=2.3, s=f"{time_scheme.__name__}:\n T={t[-1]},\nN={len(t)},\ndelta_t={T/N:.2e}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.axis("equal")
    plt.show()
    return


#plot_orbit(oscilador, Puntos_graficar, t, Euler)
#plot_orbit(oscilador, Puntos_graficar, t, RK4)
#plot_orbit(oscilador, Puntos_graficar, t, Leapfrog)
#plot_orbit(oscilador, Puntos_graficar, t, Evil_Euler, tol_N=newton_tol, tol_J=tol_jacobian, max_iter=N_max)
#plot_orbit(oscilador, Puntos_graficar, t, Crank_Nicolson, tol_N=newton_tol, tol_J=tol_jacobian, max_iter=N_max)





def stability_region(temporal_scheme):

    def R(omega, temporal_scheme):

        if temporal_scheme == Euler:
            return 1+omega  
            
        elif temporal_scheme == Crank_Nicolson:
            return (1+0.5*omega)/(1-0.5*omega)
            
        elif temporal_scheme == Evil_Euler:
            return  1/(1-omega)
        
        elif temporal_scheme == RK4:
            return 1+omega+0.5*omega**2+(1/6)*omega**3+(1/24)*omega**4
        
        elif temporal_scheme == Leapfrog:     
            return [omega + sqrt(omega**2 + 1), omega - sqrt(omega**2 + 1)]

        else:
            print("No te inventes nombres de esquemas temporales")
            return
        

    N = 2000
    real_part = uniform(-5, 5, N)
    imag_part = uniform(-5, 5, N)

    if temporal_scheme == Leapfrog:
        real_part = 0
        imag_part = uniform(-5, 5, N)

    Re, Im = meshgrid(real_part, imag_part)

    omega = Re + 1j*Im  # poner1

    if temporal_scheme != Leapfrog:
        plt.scatter(real(omega[abs(R(omega, temporal_scheme))<1]), imag(omega[abs(R(omega, temporal_scheme))<1]), s=1)  # se pintan las regiones
        
    else:
        aux = (abs(R(omega, temporal_scheme)[0]) < 1) & (abs(R(omega, temporal_scheme)[1]) < 1)
        plt.scatter(real(omega[aux]), imag(omega[aux]), s=1)

    
    plt.axis("equal")
    plt.axhline(0, color='black', linewidth=1)  
    plt.axvline(0, color='black', linewidth=1)  
    plt.xlabel(r"Re($\omega$)")
    plt.ylabel(r"Im($\omega$)")
    plt.show()


stability_region(RK4)
stability_region(Evil_Euler)
stability_region(Leapfrog)
stability_region(Crank_Nicolson)