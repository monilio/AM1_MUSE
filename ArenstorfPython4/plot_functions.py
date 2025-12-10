import matplotlib.pyplot as plt
from numpy import abs,max,sqrt
from scipy.optimize import fsolve




def plot_arenstorf(U_sol, mu):
    

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

    x_3 = U_sol[:, 0]
    y_3 = U_sol[:, 1]


    plt.figure(figsize=(8, 8))

    # 1. Establecer el color de fondo de los ejes a negro
    plt.gca().set_facecolor('black') 

    # Variables de color para visibilidad
    TEXT_COLOR = 'white'

    # Puntos y líneas originales (ajusté algunos colores para el contraste)
    plt.plot(x_3, y_3, label=' Trajectory Body 3', color='yellow') # Trayectoria en amarillo
    plt.scatter(-mu,0,          color='blue',                   s=200,  zorder=5, label='Body 1 (Earth)')
    plt.scatter(1-mu,0,         color='white',                  s=100,  zorder=5, label='Body 2 (Moon)') # Luna en blanco
    plt.scatter(L1_x, 0,        color='lime',       marker='x', s=50,   zorder=6, label='$L_1$')
    plt.scatter(L2_x, 0,        color='lime',       marker='x', s=50,   zorder=6, label='$L_2$')
    plt.scatter(L3_x, 0,        color='lime',       marker='x', s=50,   zorder=6, label='$L_3$')
    plt.scatter(L4_x, L4_y,     color='magenta',    marker='x', s=50,   zorder=6, label='$L_4$')
    plt.scatter(L5_x, L5_y,     color='magenta',    marker='x', s=50,   zorder=6, label='$L_5$')
    plt.scatter(x_3[0], y_3[0], color='red',                    s=5,    zorder=5, label='Initial position Body 3')

    # 2. Cambiar el color del texto del título y las etiquetas a blanco
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Arenstorf Orbit CR3BP')

    # 3. Cambiar el color de la rejilla a un color claro (gris)
    plt.grid(True, color='gray', linestyle='--', alpha=0.5) 

    # 4. Ajustar el color de las etiquetas de las marcas de los ejes a blanco
    plt.tick_params(axis='x')
    plt.tick_params(axis='y')

    # 5. Ajustar la leyenda para un fondo oscuro
    plt.legend(facecolor='darkgray', edgecolor='white', labelcolor=TEXT_COLOR) 

    max_abs_coord = max([abs(L3_x), abs(L4_y), abs(L5_y)])
    plot_limit = max_abs_coord * 1.1 
    plt.xlim(-plot_limit, plot_limit) 
    plt.ylim(-plot_limit, plot_limit)
    plt.axis('equal')
    plt.show() 





    plt.figure(figsize=(8, 8))

    # 1. Establecer el color de fondo de los ejes a negro
    plt.gca().set_facecolor('black') 

    # Variables de color para visibilidad
    TEXT_COLOR = 'white'

    # Puntos y líneas originales (ajusté algunos colores para el contraste)
    plt.plot(x_3, y_3, label=' Trajectory Body 3', color='yellow')              # Trayectoria en amarillo
    plt.scatter(1-mu,0, color='white', s=100, zorder=5, label='Body 2 (Moon)')  # Luna en blanco
    plt.scatter(x_3[0], y_3[0], color='red', s=10, zorder=5, label='Initial position Body 3')

    # 2. Cambiar el color del texto del título y las etiquetas a blanco
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Arenstorf Orbit CR3BP Zoom')

    # 3. Cambiar el color de la rejilla a un color claro (gris)
    plt.grid(True, color='gray', linestyle='--', alpha=0.5) 

    # 4. Ajustar el color de las etiquetas de las marcas de los ejes a blanco
    plt.tick_params(axis='x')
    plt.tick_params(axis='y')

    # 5. Ajustar la leyenda para un fondo oscuro
    plt.legend(facecolor='darkgray', edgecolor='white', labelcolor=TEXT_COLOR) 
    plt.axis('equal')
    x_initial = x_3[0]
    y_initial = y_3[0]
    zoom_half_width = 0.1
    plt.xlim(x_initial - zoom_half_width, x_initial + zoom_half_width) 
    plt.ylim(y_initial - zoom_half_width, y_initial + zoom_half_width)

    plt.show() 









"""
----------------------------
LAGRANGE POINTS FOR CIRCULAR RESTRICTED 3 BODY PROBLEM
----------------------------
d2x/dt=dx/dt=0
d2x/dt=dx/dt=0

Then: x -(1-mu)*(x+mu)/sqrt((x+mu)^2+y^2)^3 -mu*(x-1+mu)/sqrt((x-1+mu)^2+y^2)^3=0
      y -(1-mu)*y/sqrt((x+mu)^2+y^2)^3 -mu*y/sqrt((x-1+mu)^2+y^2)^3=0
"""
def Lagrange_Points(mu):
    
    def Equation(x):
        return x - (1 - mu) * (x + mu) / abs((x + mu))**3 - mu * (x - 1 + mu) / abs((x - 1 + mu))**3
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
