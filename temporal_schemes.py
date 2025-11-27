from numpy import concatenate, zeros
from aux_functions import  Newton




def Euler(F, U, t1, t2, **kwargs):

    return U + (t2-t1)*F(U, t1)





def Crank_Nicolson(F, U, t1, t2, tol_J=None, max_iter=None, tol_N=None, **kwargs):

    def crank_nicolson_function(U_next, F, U_n, t1, t2):                # funcion a resolver del Crank-Nicolson

        dt = t2-t1
        return (0.5*(F(U_next, t2) + F(U_n, t1))) - ((U_next - U_n)/dt) # 0.5*(F(U^(n+1),t(n+1)) + F(U^(n),t(n))) - (U^(n+1) - U^(n))/dt = 0

    sol = Newton(crank_nicolson_function, U, tol_J, max_iter, tol_N, F, U, t1, t2) # Se resuelve con un Newton y palante

    return sol





def RK4(F, U, t1, t2, **kwargs):
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





def Leapfrog(F, U, t1, t2):

    dt = t2-t1
    N = len(U)
    
    a = F(U, t1)[N//2:N]                                        # Se toma solo la parte de las aceleraciones

    v_mid = U[N//2:N] + 0.5*a*dt                                # v(n+1/2) = v(n) + (1/2)*a(n)*dt

    x_next = U[0:N//2] + v_mid*dt                               # x(n+1) = x(n) + v(n+1/2)*dt
    a_next = F(concatenate((x_next, zeros(N//2))), t2)[N//2:N]  # a(n+1) = F(x(n+1),t(n+1)) cogiendo solo las aceleraciones
    v_next = v_mid + 0.5*a_next*dt                              # v(n+1) = v(n+1/2) + (1/2)*a(n+1)*dt
 
    return concatenate((x_next, v_next))


