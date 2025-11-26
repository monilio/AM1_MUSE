from numpy import concatenate, array, zeros
from numpy.linalg import norm

def F(U, t):

    r = U[0:2]
    rd = U[2:4]

    return concatenate((rd, -r/norm(r)**3), axis=None) # x'' = -x/|x|^3



def oscilador(U, t):

    return array([U[1], -U[0]]) # x'' = -x



def N_body_problem(U, t, N=5):

    # dx/dt = v
    # dv/dt = sum_{j!=i}(   (x_j - x_i)/|x_j - x_i|^3   )

    # U = (pos0, vel0, pos1, vel1, ..., posN-1, velN-1)
    # Siendo posi = (x_i, y_i, z_i) y veli = (vx_i, vy_i, vz_i)


    dim = 3 # 3 componentes en cada pos y en cada vel

    F = zeros(2*N*dim)
    
    for i in range(0, N):
        
        F[(2*i)*dim:(2*i+1)*dim] = U[(2*i+1)*dim:(2*i+2)*dim]
        # F[0:3] = U[3:6]  --> dx0/dt = v0
        # F[6:9] = U[9:12] --> dx1/dt = v1
        # ...
            
        sum = zeros(dim)
        for j in range(0, N):
            if j != i:
                xj = U[2*j*dim:(2*j+1)*dim]
                xi = U[2*i*dim:(2*i+1)*dim]
                sum = sum + (xj-xi)/norm(xj-xi)**3
                
            #print(j)
            #print(i)
        F[(2*i+1)*dim:(2*i+2)*dim] = sum[:]
    
    return F
