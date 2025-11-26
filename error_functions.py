from numpy import zeros, linspace, log
from numpy.linalg import norm
from temporal_schemes import Euler, Crank_Nicolson, RK4, Evil_Euler
from aux_functions import Cauchy_problem



def richardson_extrapolation(F, U0, t1, temporal_scheme, one_of_the_cool_ones = True, **kwargs):

    N1 = len(t1)
    t2 = linspace(t1[0], t1[N1-1], 2*(N1)-1)    # Malla con el doble de puntos menos 1 y paso dt/2

    # Sacas la solucion para ambas mallas
    U_1 = Cauchy_problem(F, U0, t1, temporal_scheme, **kwargs)
    U_2 = Cauchy_problem(F, U0, t2, temporal_scheme, **kwargs)

    Err = zeros(N1)

    if one_of_the_cool_ones == True:

        if temporal_scheme == Euler:
            q = 1

        elif temporal_scheme == Crank_Nicolson:
            q = 2

        elif temporal_scheme == RK4:
            q = 4

        elif temporal_scheme == Evil_Euler:
            q = 1

        else:
            print("No te inventes nombres")
            return -1
        
        for n in range(0, N1):
            Err[n] = norm(U_1[n,:]-U_2[2*n,:])/(1-1/(2**q))

    else:
        
        for n in range(0, N1):
            Err[n] = norm(U_1[n,:]-U_2[2*n,:])
    
    return Err




def convergence_rate(F, U0, t1, temporal_scheme, N0, num_mejoras=5, **kwargs):

    N = N0

    errors = zeros(num_mejoras)
    log_N = zeros(num_mejoras)

    t = linspace(t1[0], t1[-1], N)
    U1 = Cauchy_problem(F, U0, t, temporal_scheme, **kwargs)

    for i in range(num_mejoras):

        N_fina = 2*N - 1

        t_fina = linspace(t1[0], t1[-1], N_fina)

        U2 = Cauchy_problem(F, U0, t_fina, temporal_scheme, **kwargs)

        errors[i] = norm(U2[-1,:] - U1[-1,:])
        log_N[i] = log(N_fina)

        U1 = U2

        N = N_fina

    log_errors = log(errors)

    return log_N, log_errors