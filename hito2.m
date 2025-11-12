clc, clearvarse, close all






function evil_eul = Euler_inverso(U, t, dt, F, iter)
    U_next = Euler(U, t, dt, F);
    for i = 1:iter
        U_next = U + dt*F(U_next, t+dt);
        t = t + dt;
    end

    evil_eul = U_next;

end



function runge = RK4(U, t, dt, F, iter)
    for i = 1:iter
        k1 = F(U,t);
        k2 = F(U + dt/2 * k1, t + dt/2);
        k3 = F(U + dt/2 * k2, t + dt/2);
        k4 = F(U + dt * k3, t + dt);
        U = U + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
    end

    runge = U;
end



function crank = CrankNicolson(U, t, dt, F, iter)

    U_next = Euler(U, dt, F);
    for i = 1:iter
        U_next = U + 0.5 * dt * (F(U,t) + F(U_next, t + dt));
        t = t + dt;
    end

    crank = U_next;

end


function eul = Euler(U, t, dt, F)
    eul = U + dt*F(U,t);
end