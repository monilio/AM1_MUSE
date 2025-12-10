from numpy import array,concatenate,zeros,sqrt,log, linalg, array,log10
import matplotlib.pyplot as plt


"""
----------------------------
ESTIMACIÓN DEL ERROR USANDO LA FUNCIÓN DE EXTRAPOLACIÓN DE RICHARDSON
----------------------------

Extrapolación de Richardson: phi = phi(h) + c1*h^p + c2*h^p... donde phi es la solución exacta, phi(h) es una aproximación usando un paso de tamaño h, y c1,c2... son ctes desconocidas

La idea es eliminar los términos con ctes usando dos pasos de tiempo diferentes: h y h/r

Finalmente se llega a: phi=(phi(h)-r^p*phi(h/r))/(1-r^p)

El ERROR con un tamaño de paso h/r se define como E(h/r) = phi - phi(h/r)

Introduciendo la extrapolación de Richardson: E(h/r) = (phi(h)-r^p*phi(h/r))/(1-r^p) - phi(h/r)

El orden de p para los diferentes esquemas:
- Euler/Euler Inverso: p=1
- Cranck-Nicolson: p=2
- RK4: p=4
"""
def Get_Richardson_Error(Temporal_Scheme, F, U0, tspan, N_for_h_r, p):
    
    r = 2 #División del paso
    
    _, phi_h_r_full = Temporal_Scheme(F, tspan, U0, N_for_h_r) # 1. Simulación precisa (paso h/r) -> N pasos
    phi_h_r = phi_h_r_full[-1, :] # Mantiene los últimos valores
    
    N_for_h = N_for_h_r // r
    _, phi_h_full = Temporal_Scheme(F, tspan, U0, N_for_h) # 2. Simulación menos precisa (paso h) -> N/r pasos
    phi_h = phi_h_full[-1, :] # Mantiene los últimos valores
    
    # Error de Richardson: E(h/r) = (phi(h)-r^p*phi(h/r))/(1-r^p) - phi(h/r) = (phi(h/r)-phi(h))/(r^p-1)
    error_vec = (phi_h_r - phi_h) / (r**p - 1)
    
    return linalg.norm(error_vec)



def Graficar_Convergencia_Richardson(Temporal_Scheme, F, U0, Tf_fijo, lista_N_base, p_teorico):
    
    pasos_h = []
    errores = []

    tspan = [0, Tf_fijo]
    for N in lista_N_base:
        h = Tf_fijo / N #Calcula para varios pasos N
        
        err = Get_Richardson_Error(Temporal_Scheme, F, U0, tspan, N, p_teorico) # Llama a la función de Error de Richardson
        
        pasos_h.append(h) # Almacena el nuevo h en cada paso
        errores.append(err) # Almacena el nuevo error en cada paso
        print(f"N={N:4d} | h={h:.5f} | Error Est={err:.5e}")

    #-------------------------------GRÁFICA---------------------------------------#
    plt.figure(figsize=(10, 6))
    plt.loglog(pasos_h, errores, 'bo-', lw=2, label=f'Error {Temporal_Scheme.__name__}')
    h_arr = array(pasos_h)
    C = errores[-1] / (pasos_h[-1]**p_teorico) 
    plt.loglog(h_arr, C * (h_arr**p_teorico), 'r--', label=f'Referencia Orden {p_teorico}')
    plt.grid(True, which="both", alpha=0.4)
    plt.xlabel('Paso de tiempo h (log)')
    plt.ylabel('Error Estimado Richardson (log)')
    plt.title(f'Convergencia: {Temporal_Scheme.__name__}')
    plt.legend()
    plt.gca().invert_xaxis()
    plt.show()



