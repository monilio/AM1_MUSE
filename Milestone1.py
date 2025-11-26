from numpy import array, concatenate, zeros
from numpy.linalg import norm, solve
import matplotlib.pyplot as plt

def F(U):
    
    r = U[0:2]
    v = U[2:4]
    return concatenate((v, -r/norm(r)**3), axis=None)

T = 10   
N = 10000 
dt = T/N 


######  EULER EXPLICITO  ######

U_Euler = zeros((N+1, 4))
U_Euler[0, :] = array([1, 0, 0, 1]) 

for n in range(0, N):
    U_Euler[n+1,:] = U_Euler[n,:] + dt * F(U_Euler[n, :])

plt.axis("equal")
plt.plot(U_Euler[:, 0], U_Euler[:, 1], label=f"N={N}, dt={dt}")
plt.xlabel("x")
plt.ylabel("y")
plt.show()






######  CRANK-NICOLSON  ######

U_Crank = zeros((N+1, 4))
U_Crank[0, :] = array([1, 0, 0, 1]) 

h = 1e-11

tol = 1e-11
max_iter = 1000

for n in range(0, N):

    x = U_Crank[n, :]

    for k in range(0, max_iter):
                
        J = zeros((len(U_Crank[0,:]), len(U_Crank[0,:])))
        h_vec = zeros(len(U_Crank[0,:]))
        
        for j in range(0, len(U_Crank[0,:])):
            h_vec[j] = h
            f_plus  =   0.5*(F(x+h_vec) + F(U_Crank[n,:]))-((x+h_vec-U_Crank[n,:])/dt)
            f_minus =   0.5*(F(x-h_vec) + F(U_Crank[n,:]))-((x-h_vec-U_Crank[n,:])/dt)

            J[:, j] =   (f_plus-f_minus)/(2*h)

            h_vec[j] = 0
        
        minus_f_n = -(0.5*(F(x)+F(U_Crank[n,:])) -((x-U_Crank[n,:])/dt))
        newton_diff = solve(J, minus_f_n)

        x_next = x + newton_diff
        x = x_next
        if(norm(newton_diff)<tol):
            U_Crank[n+1,:] = x
            break

plt.axis("equal")
plt.plot(U_Crank[:, 0], U_Crank[:, 1], label=f"N={N}, dt={dt}")
plt.xlabel("x")
plt.ylabel("y")
plt.show()





###### RUNGE KUTTA 4 ######

U_RK4 = zeros((N+1, 4))

U_RK4[0, :] = array([1, 0, 0, 1]) 

for n in range(0, N):

    k1 = F(U_RK4[n,:])
    k2 = F(U_RK4[n,:]+0.5*k1*dt)
    k3 = F(U_RK4[n,:]+0.5*k2*dt)
    k4 = F(U_RK4[n,:]+k3*dt)
    U_RK4[n+1,:] = U_RK4[n,:] + (1.0/6.0)*dt*(k1 + 2*(k2+k3) + k4)

plt.axis("equal")
plt.plot(U_RK4[:, 0], U_RK4[:, 1], label=f"N={N}, dt={dt}")
plt.xlabel("x")
plt.ylabel("y")
plt.show()









