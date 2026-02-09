% function for the RHS 

function rhs = rhs_advection_1d(t,U,geom,params)

  Nf  = geom.Nf;
  Nx  = geom.Nx;
  dx  = geom.dx; 
  a   = params.a;

% Initalization of Vectors    
LLF = zeros(Nf,1);
rhs = zeros(Nx,1); 

% Fill physical flux values
pf = a * U; 

% Loop for  Local Lax–Friedrichs (LLF) flux
for i = 1:Nf-1
    if i == 1 
        LLF(i) = 0.5 * (pf(Nx) + pf(i)) - 0.5 * a * (U(i)-U(Nx));
        LLF(Nf) = LLF(i); 

    else 
        LLF(i) = 0.5 * (pf(i-1) + pf(i)) - 0.5 * a * (U(i)-U(i-1));
    
    end 
    % Periodic cell average values-- Ucell(1) = Ucell(Nx) and Ucell(Nx+1) = Ucell(2)  

end 

% Loop for the residual 
for i = 1:Nx
    rhs(i) = (-1/dx(i)) * (LLF(i+1)-LLF(i));

end 
