% 1D linear advection with first order finite volume


CFL        = 0.5; 
params.CFL = CFL;
a          = 1;
params.a   = a; 
L          = 1;
params.L   = L;
Step       = [1,5,10,25,50];
T          = 0.5;
err_SU     = zeros(9,1);  % smooth uniform

for C = 1:2 

        if C == 1 % smooth case
            for J = 1:2 
                if J ==1  % Smooth Uniform case
                    for i = 1:9
                  
                        Nx = 2^(i+3);         
                        xf = linspace(0,L,Nx+1); % face locations
                       
                        
                       % Struts for geom
                         grid.Nf = Nx + 1; 
                         grid.xf = xf; 
                         
                         % Geom outputs
                         geom = GEOM_1D(grid);
                         xc = geom.xc;  
                         dx = geom.dx; 

                         % Inputs for RHS and SSPRK
                         dt = CFL*min(dx)/abs(a);
                         Nt = T/dt; 
                         t  = 0;
                         U  = sin(2*pi*xc);
                         U0 = U; 
                         % Time loop for SSPRK
                         for O = 1:Nt
                            % Call SSPRK2 fucntion
                            U_new = ssprk2(t, U, dt, @rhs_advection_1d, geom, params);
                            U = U_new;
                            t = t+dt; 
                         end 
                         % Exact solution
                         u_exact = sin(2*pi*mod(xc - a*T, L));
                         % Error values 
                         err_SU(i) = sqrt(sum((U - u_exact).^2 .* dx));
                          % Plot Smooth Uniform 
                          figure;
                          plot(xc, U0, 'k:', xc, U, 'r--', xc, u_exact, 'b-');
                          legend('u(x,0)', 'u(x,T)', 'u_{exact}(x,T)');
                          xlabel('x'); ylabel('u');
                          title(sprintf('Smooth Uniform, Nx = %d, T = %.2f', Nx, T));
                    end 
                   
                else % Smooth nonuniform case
                     for i = 1:9
                         Nx = 2^(i+3); 
                         s    = linspace(0,1,Nx+1);
                         beta = 2.0; % stretching strength
                         xf   = L*(exp(beta*s)-1)/(exp(beta)-1);
                     
                          
                         % Struts for geom
                         grid.Nf = Nx + 1; 
                         grid.xf = xf; 
                         
                         % Geom outputs
                         geom = GEOM_1D(grid);
                         xc = geom.xc;  
                         dx = geom.dx; 

                         % Inputs for RHS and SSPRK
                         dt = CFL*min(dx)/abs(a);
                         Nt = T/dt; 
                         t  = 0;
                         U  = sin(2*pi*xc);
                         U0 = U;
                         % Time loop for SSPRK
                         for O = 1:Nt
                            % Call SSPRK2 fucntion
                            U_new = ssprk2(t, U, dt, @rhs_advection_1d, geom, params);
                            U = U_new;
                            t = t+dt; 
                         end 

                         % Exact solution
                         u_exact = sin(2*pi*mod(xc - a*T, L));
                         % lLot Non-smooth Uniform:
                         figure;
                         plot(xc, U0, 'k:', xc, U, 'r--', xc, u_exact, 'b-');
                         legend('u(x,0)', 'u(x,T)', 'u_{exact}(x,T)');
                         xlabel('x'); ylabel('u');
                         title(sprintf('Smooth non-Uniform, Nx = %d, T = %.2f', Nx, T));
                   

                     end 
                     
                end 
            end

        else % non-smooth case  
            for J = 1:2
                for P = 1:5 
                    Nx = 100 * Step(P);
                    if J == 1 % Non-smooth Uniform case 
                        xf = linspace(0,L,Nx+1); % face locations
                        U  = zeros(Nx,1);
                        % Struts for geom
                         grid.Nf = Nx + 1; 
                         grid.xf = xf; 
                         
                         % Geom outputs
                         geom = GEOM_1D(grid);
                         xc = geom.xc;  
                         dx = geom.dx; 

                         % Inputs for RHS and SSPRK
                         dt = CFL*min(dx)/abs(a);
                         Nt = T/dt; 
                         t  = 0;
                         U(xc >= 0.25 & xc <= 0.5) = 1; 
                         U0 = U;
                         % Time loop for SSPRK
                         for O = 1:Nt
                            % Call SSPRK2 fucntion
                            U_new = ssprk2(t, U, dt, @rhs_advection_1d, geom, params);
                            U = U_new;
                            t = t+dt; 
                         end 
                         % Exact solution
                         x_shift = mod(xc - a*T, L);
                         u_exact = zeros(Nx,1);
                         u_exact(x_shift >= 0.25 & x_shift <= 0.5) = 1;

                         % Plot Non-smooth Uniform:
                         figure;
                         plot(xc, U0, 'k:', xc, U, 'r--', xc, u_exact, 'b-');
                         legend('u(x,0)', 'u(x,T)', 'u_{exact}(x,T)');
                         xlabel('x'); ylabel('u');
                         title(sprintf('Non-smooth Uniform, Nx = %d, T = %.2f', Nx, T));
                        
                      
                    else % Non-smooth non-uniform case
                         s      = linspace(0,1,Nx+1);
                         beta   = 2.0; % stretching strength
                         xf     = L*(exp(beta*s)-1)/(exp(beta)-1);
                         U      = zeros(Nx,1); 
                          
                        % Struts for geom
                         grid.Nf = Nx + 1; 
                         grid.xf = xf; 
                         
                         % Geom outputs
                         geom = GEOM_1D(grid);
                         xc = geom.xc;  
                         dx = geom.dx; 

                         % Inputs for RHS and SSPRK
                         dt = CFL*min(dx)/abs(a);
                         Nt = T/dt; 
                         t  = 0;
                         U(xc >= 0.25 & xc <= 0.5) = 1; 
                         U0 = U;
                         % Time loop for SSPRK
                         for O = 1:Nt
                            % Call SSPRK2 fucntion
                            U_new = ssprk2(t, U, dt, @rhs_advection_1d, geom, params);
                            U = U_new;
                            t = t+dt; 
                         end 
                         % Exact solution
                         x_shift = mod(xc - a*T, L);
                         u_exact = zeros(Nx,1);
                         u_exact(x_shift >= 0.25 & x_shift <= 0.5) = 1;

                         % Plot Non-smooth Nonuniform:
                         figure;
                         plot(xc, U0, 'k:', xc, U, 'r--', xc, u_exact, 'b-');
                         legend('u(x,0)', 'u(x,T)', 'u_{exact}(x,T)');
                         xlabel('x'); ylabel('u');
                         title(sprintf('Non-smooth Nonuniform, Nx = %d, T = %.2f', Nx, T));
               

                    end 
                end 
            end 



        end 
end 

    
                        
  Nx_vals = 2.^(4:12);
  fprintf('  Nx        ||e||_2\n');
  fprintf('-------------------\n');
  for i = 1:9
      fprintf('%5d    %.6e\n', Nx_vals(i), err_SU(i));
  end

 T_err = table(Nx_vals', err_SU, 'VariableNames', {'Nx', 'Error'});
  writetable(T_err, 'error_table.csv');