%GEOM_1D

function geom = GEOM_1D(grid)


Nf = grid.Nf; % Total number of nodes or faces
xf = grid.xf; % Pulls node cords


Nx = Nf-1; % Number of Cell Centers

% initialziation 
dx  = zeros(Nx,1); % for noninform cell widths
xcu = zeros(Nx,1); % cell centers-unifrom case
xc  = zeros(Nx,1); % cell cetners non uniform case

% dx for contstant grid
dxu = xf(2) - xf(1); 


% Builds dx and cell centers for uniform and non-uniform cases
for i = 1:Nx 
    
    I = i;

    xcu(i) = xf(I)+ 0.5*dxu; % cell centers for uniform case

    dx(i) = xf(I+1) -xf(I);  % non uniform cell width

    xc(i) = xf(I) + 0.5*dx(I); % non uniform cell centers

end


geom.Nf  = Nf;
geom.dxu = dxu;
geom.dx  = dx;
geom.xc  = xc;
geom.xcu = xcu;
geom.Nx  = Nx; 
geom.xf  = xf; 