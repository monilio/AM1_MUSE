from differential_equations import N_body_problem
from temporal_schemes import Euler, RK4, Crank_Nicolson, Evil_Euler, Leapfrog
from aux_functions import Cauchy_problem
from numpy import linspace, array
from matplotlib import pyplot as plt



def plot_NB(problem, U0_list, t, temporal_scheme, N=3, **kwargs):

    ax = plt.figure().add_subplot(projection='3d')
    dim = 3

    for U0 in U0_list:
        U = Cauchy_problem(problem, U0, t, temporal_scheme, **kwargs)

        for i in range(0, N):
            ax.plot(U[:,2*i*dim], U[:,2*i*dim+1], U[:,2*i*dim+2], label=f"body {i}")

    ax.set_title(f"Posiciones")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend()
    plt.show()
    return

T = 4
N_T = 500
t = linspace(0, T, N_T)
U0 = array([
    # r, v
     0.7, -1.1,  0.2,
    -0.15,  0.05, -0.20,

    -0.4,  0.9, -1.3,
    0.10, -0.25,  0.05,

    1.2,  0.3,  0.8,
    0.05,  0.15, -0.10,

    -1.0, -0.6,  1.1,
    -0.20,  0.00,  0.25,

    0.5,  1.4, -0.7,
    0.30, -0.10,  0.05,
])
U0 = [U0]
tol_J = 1e-8
max_iter = 1000
tol_N = 1e-9



plot_NB(N_body_problem, U0, t, RK4, tol_J=tol_J, tol_N=tol_N, max_iter=max_iter)
