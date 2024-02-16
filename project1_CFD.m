% Parameters
Lx = 2;     % Length of the domain in x-direction [uint of meter]
Ly = 1;     % Length of the domain in y-direction [uint of meter]
Lz = 2;     % Length of the domain in z-direction [uint of meter]
Nx = 65;       % Number of grid points in x-direction
Ny = 65;      % Number of grid points in y-direction
Nz = 65;       % Number of grid points in y-direction
dx = Lx / (Nx-1); % Grid spacing in x-direction
dy = Ly / (Ny-1); % Grid spacing in y-direction
dz = Lz / (Nz-1); % Grid spacing in z-direction
dt = 0.01;       % Time step size
T = 2;            % Total simulation time
nu = 0.002;       % viscosity
rho = 1000;        %density
tolerance = 0.001;%set tolerance for convergence poisson solver of Pressure field
maxIter = 100;    %maximum number of iterations allowed for poisson solver
u_inflow = 1;     %x-component velocity of inflow 
P_inflow =  0.5 * rho * u_inflow^2;

% Initialize variables
u = zeros(Nx, Ny, Nz);
v = zeros(Nx, Ny, Nz);
w = zeros(Nx, Ny, Nz);
p = ones(Nx, Ny, Nz);

% Main time loop
for t = 0:dt:T
    % Save previous velocity fields
    u_prev = u;
    v_prev = v;
    w_prev = w;


    % Predictor step (explicit)
    for i = 2:Nx-1
        for j = 2:Ny-1
           for k = 2:Nz-1 
              % Advection terms
              du_dx = (u_prev(i+1, j, k) - u_prev(i-1, j, k)) / (2*dx);
              du_dy = (u_prev(i, j+1, k) - u_prev(i, j-1, k)) / (2*dy);
              du_dz = (u_prev(i, j, k+1) - u_prev(i, j, k-1)) / (2*dz);
              
              dv_dx = (v_prev(i+1, j, k) - v_prev(i-1, j, k)) / (2*dx);
              dv_dy = (v_prev(i, j+1, k) - v_prev(i, j-1, k)) / (2*dy);
              dv_dz = (v_prev(i, j, k+1) - v_prev(i, j, k-1)) / (2*dz);
              
              dw_dx = (w_prev(i+1, j, k) - w_prev(i-1, j, k)) / (2*dx);
              dw_dy = (w_prev(i, j+1, k) - w_prev(i, j-1, k)) / (2*dy);
              dw_dz = (w_prev(i, j, k+1) - w_prev(i, j, k-1)) / (2*dz);
            
              % Viscous terms
              d2u_dx2 = (u_prev(i+1, j, k) - 2*u_prev(i, j, k) + u_prev(i-1, j, k)) / dx^2;
              d2u_dy2 = (u_prev(i, j+1,k) - 2*u_prev(i, j, k) + u_prev(i, j-1, k)) / dy^2;
              d2u_dz2 = (u_prev(i, j,k+1) - 2*u_prev(i, j, k) + u_prev(i, j, k-1)) / dz^2;
              
              d2v_dx2 = (v_prev(i+1, j, k) - 2*v_prev(i, j, k) + v_prev(i-1, j, k)) / dx^2;
              d2v_dy2 = (v_prev(i, j+1, k) - 2*v_prev(i, j, k) + v_prev(i, j-1, k)) / dy^2;
              d2v_dz2 = (v_prev(i, j, k+1) - 2*v_prev(i, j, k) + v_prev(i, j, k-1)) / dz^2;              

              d2w_dx2 = (w_prev(i+1, j, k) - 2*w_prev(i, j, k) + w_prev(i-1, j, k)) / dx^2;
              d2w_dy2 = (w_prev(i, j+1, k) - 2*w_prev(i, j, k) + w_prev(i, j-1, k)) / dy^2;
              d2w_dz2 = (w_prev(i, j, k+1) - 2*w_prev(i, j, k) + w_prev(i, j, k-1)) / dz^2; 
              
              % Predictor step
              u(i, j) = u_prev(i, j, k) - dt .* (u_prev(i, j, k).*du_dx + v_prev(i, j, k).*du_dy + w_prev(i, j, k).*du_dz - nu.*(d2u_dx2 + d2u_dy2+ d2u_dz2));
              v(i, j) = v_prev(i, j, k) - dt .* (u_prev(i, j, k).*dv_dx + v_prev(i, j, k).*dv_dy + w_prev(i, j, k).*dv_dz - nu.*(d2v_dx2 + d2v_dy2+ d2v_dz2));
              v(i, j) = w_prev(i, j, k) - dt .* (u_prev(i, j, k).*dw_dx + v_prev(i, j, k).*dw_dy + w_prev(i, j, k).*dw_dz - nu.*(d2w_dx2 + d2w_dy2+ d2w_dz2));
           end
        end
    end
    
 % Apply boundary conditions 
    u(1, (Ny + 1) / 2:end,:) = u_inflow; 
    u(1, 1:(Ny + 1) / 2, :) = 0;
    u(end, 1:(Ny + 1) / 2, :) = u(end-1, 1:(Ny + 1)/ 2, :);
    u(end, (Ny + 1) / 2:end, :) =0; 
    u(1, 1:(Ny+1)/2, :) = 0;
    u(:, 1, :) = 0;   
    u(:, end, :) = 0;   
    u(:, :, 1) = 0; 
    u(:, :, end) = 0; 

    v(1, 1:(Ny + 1) / 2, :) =  -v(2, 1:(Ny + 1) / 2, :);
    v(1, (Ny + 1) / 2:end,:) = 0;
    v(end, 1:(Ny + 1) / 2 , :) = v(end-1, 1:(Ny + 1)/ 2, :); 
    v(end, (Ny + 1) / 2:end, :) =0;
    v(:, 1, :) = 0;   
    v(:, end, :) = 0;
    v(:, :, 1) = 0;   
    v(:, :, end) = 0;
   
    w(1, 1:(Ny + 1) / 2, :) = -w(2, 1:(Ny + 1) / 2, :);
    w(1, (Ny + 1) / 2:end,:) = 0;   
    w(end, 1:(Ny + 1) / 2 , :) = w(end-1, 1:(Ny + 1)/ 2, :);
    w(end, (Ny + 1) / 2:end, :) =0;
    w(:, 1, :) = 0;   
    w(:, end, :) = 0;
    w(:, :, 1) = 0;   
    w(:, :, end) = 0;
    
 
    %{ Zero Gradient (Neumann) pressure condition at outflow boundaries
    p(end, (Ny + 1) / 2:end,:) = p(end-1, (Ny + 1) / 2:end,:);  % Right outflow boundary
    p(end, 1:(Ny + 1) / 2, :) =0;
    p(1, (Ny + 1) / 2:end,:) = (0.9999)*p(2, (Ny + 1) / 2:end,:) + (0.0001) * P_inflow;   % Left inflow boundary (P_inflow is a constant)
    p(1, 1:(Ny + 1) / 2, :) = p(2, 1:(Ny + 1) / 2, :);
    p(:, end, :) = p(:, end-1, :);
    p(:, 1, :) = p(:, 2, :);
    p(:, :, end) = p(:, :, end-1);
    p(:, :, 1) = p(:, :, 2);
    % Corrector step (implicit pressure correction)
    % Use a simple Gauss-Seidel solver for simplicity
    % Iterations for pressure correction
    % Solve pressure Poisson equation
  for iter = 1:maxIter
        p_prev = p;    
        for i = 2:Nx-1
            for j = 2:Ny-1
                for k = 2:Nz-1

                    du_dx = (u(i+1, j, k) - u(i-1, j, k)) / (2*dx);
                    dv_dy = (v(i, j+1, k) - v(i, j-1, k)) / (2*dy);
                    dw_dz = (w(i, j, k+1) - w(i, j, k-1)) / (2*dz);

                    rhs_ijk = -rho .* (du_dx + dv_dy + dw_dz);

                    p(i, j, k) = (dy^2 .* (p(i+1, j, k) + p(i-1, j, k)) + dx^2 .* (p(i, j+1, k) + p(i, j-1, k)) + dz^2 * (p(i, j, k+1) + p(i, j, k-1)) + dx^2 * dy^2 * dz^2 * rhs_ijk) / (2*(dx^2 + dy^2 + dz^2));
                end
            end
        end
        
        % Check convergence
        if max(max(abs(p_prev - p))) < tolerance
            break;
        end
    end

    p_new = p;


    % Update velocities with pressure correction
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 2:Nz-1
                % Correct u-velocity
                u(i, j, k) = u(i, j, k) - dt / rho * (p_new(i+1, j, k) - p_new(i-1, j, k)) / (2 * dx);
                
                % Correct v-velocity
                v(i, j, k) = v(i, j, k) - dt / rho * (p_new(i, j+1, k) - p_new(i, j-1, k)) / (2 * dy);
                
                % Correct w-velocity
                w(i, j, k) = w(i, j, k) - dt / rho * (p_new(i, j, k+1) - p_new(i, j, k-1)) / (2 * dz);
            end
        end
    end
 

    % Apply boundary conditions 
    u(1, (Ny + 1) / 2:end,:) = u_inflow; 
    u(1, 1:(Ny + 1) / 2, :) = 0;
    u(end, 1:(Ny + 1) / 2, :) = u(end-1, 1:(Ny + 1)/ 2, :);
    u(end, (Ny + 1) / 2:end, :) =0; 
    u(1, 1:(Ny+1)/2, :) = 0;
    u(:, 1, :) = 0;   
    u(:, end, :) = 0;   
    u(:, :, 1) = 0; 
    u(:, :, end) = 0; 

    v(1, 1:(Ny + 1) / 2, :) =  -v(2, 1:(Ny + 1) / 2, :);
    v(1, (Ny + 1) / 2:end,:) = 0;
    v(end, 1:(Ny + 1) / 2 , :) = v(end-1, 1:(Ny + 1)/ 2, :); 
    v(end, (Ny + 1) / 2:end, :) =0;
    v(:, 1, :) = 0;   
    v(:, end, :) = 0;
    v(:, :, 1) = 0;   
    v(:, :, end) = 0;
   
    w(1, 1:(Ny + 1) / 2, :) = -w(2, 1:(Ny + 1) / 2, :);
    w(1, (Ny + 1) / 2:end,:) = 0;   
    w(end, 1:(Ny + 1) / 2 , :) = w(end-1, 1:(Ny + 1)/ 2, :);
    w(end, (Ny + 1) / 2:end, :) =0;
    w(:, 1, :) = 0;   
    w(:, end, :) = 0;
    w(:, :, 1) = 0;   
    w(:, :, end) = 0;
    
 
    %{ Zero Gradient (Neumann) pressure condition at outflow boundaries
    p(end, (Ny + 1) / 2:end,:) = p(end-1, (Ny + 1) / 2:end,:);  % Right outflow boundary
    p(end, 1:(Ny + 1) / 2, :) =0;
    p(1, (Ny + 1) / 2:end,:) = (0.9999)*p(2, (Ny + 1) / 2:end,:) + (0.0001) * P_inflow;;    % Left inflow boundary (P_inflow is a constant)
    p(1, 1:(Ny + 1) / 2, :) = p(2, 1:(Ny + 1) / 2, :);
    p(:, end, :) = p(:, end-1, :);
    p(:, 1, :) = p(:, 2, :);
    p(:, :, end) = p(:, :, end-1);
    p(:, :, 1) = p(:, :, 2);
    
t    
end
x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny)';
z=ones(1,1,Nz);
z(1,1,:)=linspace(0,Lz,Nz);
ux = permute(u, [2, 1, 3]);
uy = permute(v, [2, 1, 3]);
uz = permute(w, [2, 1, 3]);
N = 1000;  xstart = max(x)*rand(N,1);  ystart = max(y)*rand(N,1);  zstart = max(z)*rand(N,1);
figure(1) ;  h=streamline(x,y,z,ux,uy,uz,xstart,ystart,zstart);
title('streamline 2D ');
figure(2) ;  h=streamline(x,y,z,ux,uy,uz,xstart,ystart,zstart);
colorbar('peer',gca,'SouthOutside');
axis image;
view(3);
title('streamline 3D ');

figure(3)
contourf(x,y,uy(:,:,(Nz-1)/2));
axis image;
colorbar('peer',gca,'SouthOutside'); 
title('velosity of y_direction (m/s)');
xlabel('X [m/s]');
ylabel('Y [m/s]');

figure(4)
contourf(x,y,uz(:,:,(Nz-1)/2));
axis image;
colorbar('peer',gca,'SouthOutside'); 
title('velosity of z_direction(m/s)');
xlabel('X [m/s]');
ylabel('Y [m/s]');

figure(5)
a= sqrt(ux.^2 + uy.^2 + uz.^2);
plot(x,a(:,(Ny-1)/2,(Nz-1)/2));
title('velosity (m/s) in middle of y & Z'); hold on
plot(x,u(:,(Ny-1)/2,(Nz-1)/2)); hold on
plot(x,v(:,(Ny-1)/2,(Nz-1)/2)); 
xlabel('X [m/s]');
ylabel('velosity');
legend('velocity','x-velocity', 'y-velocity', 'Location', 'best');



figure(6)
plot(a((Nx-1)/2,:,(Nz-1)/2),y); 
title('velosity (m/s) in middle of x & Z'); hold on
plot(u((Nx-1)/2,:,(Nz-1)/2),y); hold on
plot(v((Nx-1)/2,:,(Nz-1)/2),y); hold on
xlabel('velosity');
ylabel(' Y [m/s]');
legend('velocity','x-velocity', 'y-velocity', 'Location', 'best');


figure(7)
plot(x,p(:,(Ny-1)/2,(Nz-1)/2));
title('pressure (pa) in middle of y & z axis');
xlabel('X [m]');
ylabel('pressure ');

figure(8)
plot(p((Nx-1)/2,:,(Nz-1)/2),y);
title('pressure (pa) in middle of x & z axis');
xlabel('pressure');
ylabel('Y [m]');


figure(9)
contourf(x,y,p(:,:,(Nz-1)/2));
axis image;
colorbar('peer',gca,'SouthOutside'); 
title('pressure (pa)');
xlabel('X [m]');
ylabel('Y [m]');

