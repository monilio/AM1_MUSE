from numpy import array,  zeros
from numpy.linalg import norm, solve, LinAlgError


def Jacobian(f, x, h, F, U, t1, t2):

    n = len(x)
    
    m = len(f(x, F, U, t1, t2)) # Longitud del vector que devuelve f

    J = zeros((m, n))

    dx = zeros(n)

    for i in range(0, n):
        dx[i] = h
        J[:, i] = (f(x+dx, F, U, t1, t2) - f(x-dx, F, U, t1, t2))/(2*h)
        dx[i] = 0

    return J
    




def Newton(f, x, h, N, newton_tol, F, U, t1, t2):

    x_n = array(x, copy=True)
    for n in range(0, N):
        J_n = Jacobian(f, x_n, h, F, U, t1, t2)
        f_n = f(x_n, F, U, t1, t2)
        try:
            delta_x = solve(J_n, -f_n)
        except LinAlgError:
            print("Big bad, Newton big sad")
            break
        
        x_n = x_n + delta_x
        
        if(norm(delta_x)<newton_tol):
            break

    return x_n





def Cauchy_problem(F, U0, t, temporal_scheme, **kwargs):

    k = 0
    N = len(t)
    N2 = len(U0)
    U = zeros((N, N2))
    U[0,:] = U0

    for n in range(0, N-1):
        U[n+1,:] = temporal_scheme(F, U[n,:], t[n], t[n+1], **kwargs)
        print(k)
        k = k+1

    return U