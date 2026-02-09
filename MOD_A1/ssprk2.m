% SSPRk2 function 

function U_new = ssprk2(t, U, dt, rhs, geom, params)


R = rhs(t,U,geom,params);

FE    = U + dt * R; 

U_1   = FE;

R2    = rhs(t + dt, U_1, geom, params);

FE2   = U_1 + dt*R2;

U_new = 0.5 * U + 0.5 * FE2; 
