%GEOM_1D

function geom = GEOM_1D(grid)


Nf = grid.Nf; % Total number of nodes or faces
xf = grid.xf; % Pulls node cords

Nx = Nf-1; % Number of Cell Centers

% initialziation 
dx  = zeros(Nx,1); % for noninform cell widths
xcu = zeros(Nx,1); % cell centers-unifrom case
xc  = zeros(Nx,1); % cell cetners non uniform case



% Builds dx and cell centers for uniform and non-uniform cases
for i = 1:Nx 
    
    I = i;
    dx(i)  = xf(I+1) -xf(I);  % non uniform cell width
    
    xcu(i) = xf(I)+ 0.5*dx(i); % cell centers for uniform case

    xc(i)  = xf(I) + 0.5*dx(I); % non uniform cell centers

end


geom.Nf  = Nf;
geom.dx  = dx;
geom.xc  = xc;
geom.xcu = xcu;
geom.Nx  = Nx; 
geom.xf  = xf; 