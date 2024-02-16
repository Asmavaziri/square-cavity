

% Generate cartesian grid using the meshgrid function
    Nx = 129;                        % number of point in x_direction
    Ny = 129;                        % number of point in x_direction
    Nz = 1;                         % number of point in x_direction
    theta = 35;                      % angle of ramp 
    [x y] = gridgen(Nx,Ny,theta);
    z = zeros(Nx,Ny);
    disp('Grid generated');

% grid metrics

z_width = [0,0,1];                   % Unit vector in z-dir
Ncx = Nx-1;                         % Number of cells (pts-1) in i dir
NCy = Ny-1;                         % Number of cells (pts-1) in j dir

for i = 1:Ncx
    for j = 1:NCy
        % Assemble the face lengths
        lenx_R = x(i+1,j+1)-x(i+1,j);
        leny_R = y(i+1,j+1)-y(i+1,j);
        
        lenx_T = x(i,j+1)-x(i+1,j+1);
        leny_T = y(i,j+1)-y(i+1,j+1);

        lenx_L = x(i,j)-x(i,j+1);
        leny_L = y(i,j)-y(i,j+1);

        lenx_B = x(i+1,j)-x(i,j);
        leny_B = y(i+1,j)-y(i,j);
        
        % Compute midpoint of cell (for plotting)
        xmid(i,j) = (x(i,j) + x(i+1,j))/2;
        ymid(i,j) = (y(i,j) + y(i,j+1))/2;
        
        % Compute volume of cell using A.BxC (volume of parallelepiped)
        volume(i,j) = abs(dot(z_width,cross(-1*[lenx_B,leny_B,0],...
                      [lenx_R,leny_R,0])));
        
        % Compute area of cell
        sR(i,j) = sqrt((lenx_R)^2 + (leny_R)^2);
        sT(i,j) = sqrt((lenx_T)^2 + (leny_T)^2);
        sL(i,j) = sqrt((lenx_L)^2 + (leny_L)^2);
        sB(i,j) = sqrt((lenx_B)^2 + (leny_B)^2);
        
        % Compute outward normal of faces (return 3 component vector)
        temp_nE = cross([lenx_R,leny_R, 0]/sR(i,j), z_width);
        temp_nN = cross([lenx_T,leny_T, 0]/sT(i,j), z_width);
        temp_nW = cross([lenx_L,leny_L, 0]/sL(i,j), z_width);
        temp_nS = cross([lenx_B,leny_B, 0]/sB(i,j), z_width);
        
        % Truncate normal vector to 2 components
        nE{i,j} = [temp_nE(1) temp_nE(2)];
        nN{i,j} = [temp_nN(1) temp_nN(2)];
        nW{i,j} = [temp_nW(1) temp_nW(2)];
        nS{i,j} = [temp_nS(1) temp_nS(2)];
        
        % Clear unecessary variables
        clear temp_nE temp_nN temp_nW temp_nS
        clear lenx_R lenx_B lenx_L lenx_T
        clear leny_R leny_B leny_L leny_T
    end 
end
disp('Grid Metrics calculated');


% freestream condition 
free_rho      = 1.2;                            % kg/m^3
free_T        = 300.00;                         % Kelvin
free_P        = 100000;                         % Pa
free_M        = 1.2;                            % Mach
free_aoa      = 0;                              % deg

free_a        = speedsound(free_P,free_rho);    % m/s
free_u        = free_M*free_a*cosd(free_aoa);   % m/s
free_v        = free_M*free_a*sind(free_aoa);   % m/s
free_vel      = [free_u free_v];

% Freestream Primitive State Vector (PSV)
V_free        = [free_rho free_u free_v free_P];

diary('output.txt')
fprintf('Freestream Conditions:\n');
fprintf('Mach:         %10.2f\n',free_M);
fprintf('u Velocity:   %10.2f m/s\n',free_u);
fprintf('v Velocity:   %10.2f m/s\n',free_v);
fprintf('Pressure:     %10.2e Pa\n',free_P);
fprintf('Temperature:  %10.2f K\n',free_T);
fprintf('Density:      %10.2f kg/m^3\n\n',free_rho);

%set itration variables
% Iteration variables
iterations   = 100;             % Number of iterations to run

gbl_timestep = 0;               % 0 = local timestepping
                                % 1 = global timestepping
                                
timestep     = 1e-5;            % Timestep for global timestepping
CFL          = 0.5;             % Courant number

m_stage      = 4;               % m-stage time stepping
                                % e.g. 1 for Euler step
                                %      4 for 4th-order RK

% Output variables
freq         = 10;              % Reporting frequency for output

plotcontours = 0;               % Plot contours during iterations
                                % 0 - off
                                % 1 - Density
                                % 2 - U Velocity
                                % 3 - V Velocity
                                % 4 - Pressure

plots_on     = 1 ;              % Plot PSV contour plots after 
                                % completion
                                

% set initial condition
tic
resid_i      = 0;               % Iterative residual
resid_0      = 0;               % Step 0 residual
start_iter   = 0;               % Used for multiple runs
end_iter     = 0;               % Used for multiple runs
residReduced = 0;               % If divergence detected
fm = 1;                         % Used for capturing movies

% Combine normals and areas into big cell array which will be passed
% to the function which computes the residual
normals = {nE nN nW nS};
areas   = {sR sT sL sB};

% Initalize variables which will allow for visualization
% i.e. Plot Contours
con_density = zeros([Ncx,NCy]);
con_uvel    = zeros([Ncx,NCy]);
con_vvel    = zeros([Ncx,NCy]);
con_pres    = zeros([Ncx,NCy]);

% Loop through all cells and init PSV to freestream conditions
% Convert PSV to conservative state vector
% Init residual to 0
for i = 1:Ncx
    for j = 1:NCy
        V{i,j}     = V_free;
        U{i,j}     = convV_U(V{i,j});
        resid{i,j} = [0 0 0 0];
    end
end
diary off
disp('Solution initialized');


% timestep loop
start_iter = start_iter + 1;        % Start at iteration 1
diary on

% Create plot window
if plotcontours == 1
    hc = figure('name','Density Contour'); 
elseif plotcontours == 2
    hc = figure('name','U Velocity Contour'); 
elseif plotcontours == 3
    hc = figure('name','V Velocity Contour'); 
elseif plotcontours == 4
    hc = figure('name','Pressure Contour'); 
end


% Main loop in time
for iter = start_iter:(end_iter + iterations)
    % Time variable used to measure time/iteration
    ti1 = cputime;
    
    % Initialize iteration residual to 0
    resid_i = 0;
    
    % Save CSV from this timestep (to be used in m-stage)
    U0 = U;
    
    % M-stage timestepping scheme
    for m = 1:m_stage
        % Calculate residual using function calcResid
        % Passes PSV, normals, areas cell array,
        %        freestream PSV, and nci and ncj
        resid = calcResid(V, V_free, normals, areas, Ncx, NCy);
        
        % Loop through all cells to update solution
        for i = 1:Ncx
            for j = 1:NCy
                
                % If local timestepping, calculate timestep of the cell
                if gbl_timestep == 0
                    vel = [V{i,j}(2) V{i,j}(3)];
                    cell_a = speedsound(V{i,j}(4),V{i,j}(1));
                    dt(1) = CFL * sR(i,j)/(abs(vel(1)*nE{i,j}(1) +...
                            vel(2)*nE{i,j}(2))+cell_a);
                    dt(2) = CFL * sT(i,j)/(abs(vel(1)*nN{i,j}(1) +...
                            vel(2)*nN{i,j}(2))+cell_a);
                    dt(3) = CFL * sL(i,j)/(abs(vel(1)*nW{i,j}(1) +...
                            vel(2)*nW{i,j}(2))+cell_a);
                    dt(4) = CFL * sB(i,j)/(abs(vel(1)*nS{i,j}(1) +...
                            vel(2)*nS{i,j}(2))+cell_a);
                    timestep = min(dt);
                end
                
                % Update solution using the saved CSV
                % Multiply by 'alpha' constant 
                U{i,j} = U0{i,j} -...
                         1/(m_stage-(m-1))*timestep/volume(i,j)*resid{i,j};
                
                % Update cell PSV
                V{i,j} = convU_V(U{i,j});
            
                % Update contour arrays used for plotting
                con_density(i,j) = V{i,j}(1);
                con_uvel(i,j)    = V{i,j}(2);
                con_vvel(i,j)    = V{i,j}(3);
                con_pres(i,j)    = V{i,j}(4);

                % Assemble first part of L2 norm for residual
                resid_i = resid_i + resid{i,j}.^2;
            end
        end
    end   
    
    % Assemble second part of L2 norm for residual
    resid_i = (resid_i).^.5/(Ncx*NCy);
    
    % Assign normalization value in first interation
    if iter == 1
        resid_0 = resid_i;
        
    end

    if isnan(resid_i/resid_0)
        break;
        disp('Solution corrupt.');
    end
    
    % Detects divergence happening in x-mom resid and cuts CFL in half
	if ((resid_i(2)/resid_0(2)) >= (1e1))
        if residReduced == 0
            CFL = CFL/2; 
            notice = sprintf('Divergence detected.  CFL reduced to %5.2f',CFL);
            disp(notice);
            residReduced = residReduced + 1;

        end
    end
    
    % Computes time/iteration
    ti2 = cputime-ti1;
    if mod(iter,freq) == 0
        % Plot contours if wanted
        if plotcontours == 1
            % Assembles contour plot
            [C,h] = contourf(xmid',ymid',con_density');
            
            % Removes lines from contour plot
            set(h, 'LineStyle','none');
            
            % Adds a title to the plot
            s = sprintf('Density Contour at %4d iterations',iter);
            title(s)
            
            % Adds gridlines, corrects aspect ratio
            grid on;axis image;drawnow;
            
            % Assembles array for movie viewing
            mov(fm) = getframe(gca);
            fm=fm+1;
        elseif plotcontours == 2
            [C,h] = contourf(xmid',ymid',con_uvel');
            set(h, 'LineStyle','none');
            s = sprintf('U-Vel. Contour at %4d iterations',iter);
            title(s)
            colorbar('peer',gca,'SouthOutside'); 
            grid on;axis image;drawnow;
            mov(fm) = getframe(gca);
            fm=fm+1;
        elseif plotcontours == 3
            [C,h] = contourf(xmid',ymid',con_vvel');
            set(h, 'LineStyle','none');
            s = sprintf('V-Vel. Contour at %4d iterations',iter);
            title(s)
            colorbar('peer',gca,'SouthOutside'); 
            grid on;axis image;drawnow;
            mov(fm) = getframe(gca);
            fm=fm+1;
        elseif plotcontours == 4
            [C,h] = contourf(xmid',ymid',con_pres');
            set(h, 'LineStyle','none');
            s = sprintf('Pressure Contour at %4d iterations',iter);
            title(s)
            colorbar('peer',gca,'SouthOutside'); 
            grid on;axis image;drawnow;
            mov(fm) = getframe(gcf);
            fm=fm+1;
        end

    end
end

start_iter = iter;
end_iter = iter;
diary off


if plots_on == 1
    
    
    figure 
    contourf(xmid',ymid',con_density');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('Density (kg/m^3)')

    figure 
    contourf(xmid',ymid',con_uvel');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('U Velocity (m/s)')

    figure 
    contourf(xmid',ymid',con_vvel');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('V Velocity (m/s)')

    figure 
    contourf(xmid',ymid',con_pres');
    axis image;
    colorbar('peer',gca,'SouthOutside'); 
    title('Pressure (Pa)')
    
    figure
    x1=x(1:end-1,1)';
    y1=y(1:end-1,:)';
    y1=y1(1:end-1,:);
    u1=con_uvel';
    v1=con_vvel';
    N = 1000;  xstart = max(x1).*rand(128,1);  ystart = y1(1,:).*rand(128,1);
    h=streamline(x1,y1,u1,v1,xstart,ystart,[0.1, 300]);
end

%%GENERATE GRID
function [X ,Y] = gridgen(npi,npj,theta);
meshxxx= linspace(0,6,npi);
meshyyy= linspace(0,(tan(deg2rad(theta))+5),npj);
x = ones(npi,npj) .* meshxxx' ;
y = ones(npi,npj).* meshyyy ;
theta1 = theta /((npj-1)/2);
theta2 = [0:theta1:theta];
[n,m] = size(theta2) ;
theta3 =[theta2 ones(1,(npj-m)).* theta];
kk=1;
  for o=1:npj
       for p = ((npi-1)/2)+2:npi
           y(p,o) = y(((npi-1)/2)+1,o) + (tan(deg2rad(theta-theta3(1,o))))*(6/(npi-1))*kk;
           kk=kk+1;
       end
       kk=1;
  end
  
X=x;
Y=y;
end



%%% CALCULATE FULL FLUX VECTOR 
function f = flux(V,n)
% Inputs: Primitive state vector, V, from a cell
%         Normal vector, n, from a cell face

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

%conV  = dot([u v],n);           % Contravariant velocity
conV  = u*n(1) + v*n(2);

% Assemble flux vector
f(1) = rho*conV;
f(2) = rho*u*conV + P*n(1);
f(3) = rho*v*conV + P*n(2);
f(4) = rho*h_0(V)*conV;
end


%CALCULATE NEGATIVE FLUX VECTOR
function f = fneg(V,n)
% Inputs: Primitive state vector, V, from a cell
%         Normal vector, n, from a cell face

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

gamma = 1.4;                    % Gamma
a     = speedsound(P,rho);      % Speed of sound
%conV  = dot([u v],n);           % Contravariant velocity
conV  = u*n(1) + v*n(2);
M     = conV/a;                 % Contravarient Mach number

% Assemble flux vector
f(1)  = -rho*a/4*(M-1)^2;
f(2) = f(1) * (u + n(1)*(-conV-2*a)/gamma);
f(3) = f(1) * (v + n(2)*(-conV-2*a)/gamma);
f(4) = f(1) * (h_0(V) - a^2*(M+1)^2/(gamma + 1));
end

%CALCULATE POSITIVE FLUX VECTOR
function f = fpos(V,n)
% Inputs: Primitive state vector, V, from a cell
%         Normal vector, n, from a cell face

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

gamma = 1.4;                    % Gamma
a     = speedsound(P,rho);      % Speed of sound
%conV  = dot([u v],n);           % Contravariant velocity
conV  = u*n(1) + v*n(2);
M     = conV/a;                 % Contravarient Mach number

% Assemble flux vector
f(1)  = rho*a/4*(M+1)^2;
f(2) = f(1) * (u + n(1)*(-conV+2*a)/gamma);
f(3) = f(1) * (v + n(2)*(-conV+2*a)/gamma);
f(4) = f(1) * (h_0(V) - a^2*(M-1)^2/(gamma + 1));
end


% CALCULATE SPEED OF SOUND
function c = speedsound(P,rho)
% Inputs: Pressure and Density

% speed of sound = sqrt(gamma*pressure/rho)
c = sqrt(1.4*P/rho);
end

% CALCULATE TOTAL ENTHALPY 
function h = h_0(V)
% Inputs: Primitive state vector from a cell

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% Total enthalpy = total energy + pressure/rho
h = e_0(V)+P/rho;
end

% CALCULATE TOTAL ENERGY 
function e = e_0(V)
% Inputs: Primitive state vector from a cell

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% Total energy = pressure/[(gamma-1)*rho] + 1/2(u^2+v^2)
e = P/((1.4-1)*rho)+.5*(u^2+v^2)
end

% CALCULATE PRIMITIVE STATE VECTOR 
function V = convU_V(U)
% Inputs: Conservative state vector from a cell

% Parse out variables and apply more names
rho   = U(1);   % Density
rho_u = U(2);   % Density * u vel.
rho_v = U(3);   % Density * v vel.
rho_e = U(4);   % Density * total energy

% V = [rho u_vel v_vel pressure]
V(1) = rho;
V(2) = rho_u/rho;
V(3) = rho_v/rho;
V(4) = (rho_e - rho/2*(V(2)^2+V(3)^2))*(1.4-1);
end


% CALCULATE CONSERVATIVE STATE VECTOR
function U = convV_U(V)
% Inputs: Primitive state vector from a cell

% Parse out variables and apply more names
rho = V(1);     % Density
u   = V(2);     % u velocity
v   = V(3);     % v velocity
P   = V(4);     % Pressure

% U = [rho rho*u_vel rho*v_vel rho*energy]
U(1) = rho;
U(2) = rho*u;
U(3) = rho*v;
U(4) = rho*e_0(V);
end


% CALCULATE RESIDUAL OF GRID 
function resid = calcResid(V, V_free, normals, areas, nci, ncj)
% Inputs: V       - Cell array containing the primitive state vector for
%                   each cell in the grid
%         V_free  - Primitive state vector containing freestream conditions
%         normals - Cell array containing all of the face normals for each
%                   cell
%         areas   - Cell array containing all of the face areas for each
%                   cell
%         nci,ncj - Number of cell in i and j

% Parse out normal vectors in smaller, specialized cell arrays
nE = normals{1};
nN = normals{2};
nW = normals{3};
nS = normals{4};

% Parse out areas vectors in smaller, specialized cell arrays
sE = areas{1};
sN = areas{2};
sW = areas{3};
sS = areas{4};

% Initialize cell array containing the flux for each face
e_flux = cell(nci,ncj);
w_flux = cell(nci,ncj);
n_flux = cell(nci,ncj);
s_flux = cell(nci,ncj);

% Initialize residual cell array
resid = cell(nci,ncj);

% Assemble freestream velocity vector
free_vel = [V_free(2) V_free(3)];

% Calculate freestream Mach number
free_a   = speedsound(V_free(4), V_free(1));

% Loop through all cells
for i = 1:nci
    for j = 1:ncj
        % Calculate speed of sound for cell
        cell_a = speedsound(V{i,j}(4),V{i,j}(1));
        
        % Assemble velocity vector of cell
        vel    = [V{i,j}(2) V{i,j}(3)];
        
        % Contravariant Mach for face
        %conM = dot(vel,nE{i,j})/cell_a;
        conM = (vel(1)*nE{i,j}(1) + vel(2)*nE{i,j}(2))/cell_a;
        
        % If east face is at the outflow boundary
        if i == nci
            % Boundary normal equal to face normal
            nB = nE{i,j};
            
            % Check to see is contravariant M >= 1
            if conM >= 1
                % Assign boundary flux to nci cell flux
                Vb = V{i,j};
                
                % Calculate east face flux
                e_flux{i,j} = flux(Vb,nB);
            else
                % Calculate positive & negative Riemann invariants
                rpos = (vel.*nB) + 2/(1.4-1)*cell_a;
                rneg = (free_vel.*nB) - 2/(1.4-1)*free_a;

                % Average invariants to get normal velocities
                un = (rpos+rneg)/2;

                % Obtain velocities at boundaries
                ub = vel + (un - vel).*nB;                 

                % Vb is primitive state vector at outflow boundary
                % Density is extrapolated
                % u,v are calculated
                % Pressure is extrapolated
                Vb = [V{i,j}(1) ub(1) ub(2) V{i,j}(4)];

                % Calculate East face flux
                e_flux{i,j} = flux(Vb,nB);  
            end
        % If east face is internal
        else
            % Check to see is contravariant M >= 1
            if conM >= 1
                % Calculate east face flux using full flux
                e_flux{i,j} = flux(V{i,j},nE{i,j});
            else
                % Calculate east face flux using f+ and f-
                e_flux{i,j} = fpos(V{i,j},nE{i,j}) + fneg(V{i+1,j},nE{i,j});
            end
        end
        
        
        % Contravariant Mach for face
        %conM = dot(vel,nN{i,j})/cell_a;
        conM = (vel(1)*nN{i,j}(1) + vel(2)*nN{i,j}(2))/cell_a;
        
        % If north face is invicid wall
        if j == ncj
            % Boundary normal equal to face normal
            nB = nN{i,j};
            
            % Contravariant velocity:
            % Dot product of the u,v velocity of cell(i,ncj)
            % and the wall normal vector
            conV = dot(vel, nB);
            
            % Wall velocity:
            % Velocity of cell(i,ncj) - contravariant vel. * wall normal
            % vector
            velB = vel - conV*nB;
            
            % Vb is primitive state vector at the wall
            % Density is extrapolated
            % u,v are calculated
            % Pressure is extrapolated
            Vb = [V{i,j}(1) velB(1) velB(2) V{i,j}(4)];
            
            % Calculate north face flux
            n_flux{i,j} = flux(Vb,nB);
        else
            % Check to see is contravariant M >= 1
            if conM >=1
                % Calculate north face flux using full flux
                n_flux{i,j} = flux(V{i,j},nN{i,j});
            else
                % Calculate north face flux using f+ and f-
                n_flux{i,j} = fpos(V{i,j},nN{i,j}) + fneg(V{i,j+1},nN{i,j});
            end
        end
        
        
        % Contravariant Mach for face
        %conM = dot(vel,nW{i,j})/cell_a;
        conM = (vel(1)*nW{i,j}(1) + vel(2)*nW{i,j}(2))/cell_a;
        
        % If west face is inflow boundary
        if i == 1
            % Boundary normal equal to face normal
            nB = nW{i,j}; 
            
            if conM >= 1
                % Assign boundary state vec to freestream vals
                Vb = V_free;

                % Calc flux for west face
                w_flux{i,j} = flux(Vb,nB);
            else
                % Calculate positive & negative Riemann invariant
                rpos = (vel.*nB) + 2/(1.4-1)*cell_a;
                rneg = (free_vel.*nB) - 2/(1.4-1)*free_a;

                % Average invariants to get normal velocities
                un = (rpos+rneg)/2;

                % Obtain velocities at boundaries
                ub = free_vel + (un + free_vel).*nB;

                % Vb is primitive state vector at inflow boundary
                % Density is freestream
                % u,v are calculated
                % Pressure is freesteam
                Vb = [V_free(1) ub(1) ub(2) V_free(4)];

                % Calculate west face flux
                w_flux{i,j} = flux(Vb,nB);
            end
        else
            % Check to see is contravariant M >= 1
            if conM >=1
                % Calculate west face flux using full flux
                w_flux{i,j} = flux(V{i,j},nW{i,j});
            else
                % Calculate west face flux using the negative
                % of the adjacent cell east face flux
                w_flux{i,j} = -1*e_flux{i-1,j};
            end
        end
        
        
        % Contravariant Mach for face
        %conM = dot(vel,nS{i,j})/cell_a;
        conM = (vel(1)*nS{i,j}(1) + vel(2)*nS{i,j}(2))/cell_a;
        
        % If south face is invicid wall
        if j == 1
            % Boundary normal equal to face normal
            nB = nS{i,j};
            
            % Contravariant velocity:
            % Dot product of the u,v velocity of cell(i,1)
            % and the wall normal vector
            conV = dot(vel, nB);
            
            % Wall velocity:
            % Velocity of cell(i,1) - contravariant vel. * wall normal
            % vector
            velB = vel - conV*nB;
            
            % Vb is primitive state vector at the wall
            % Density is extrapolated
            % u,v are calculated
            % Pressure is extrapolated
            Vb = [V{i,j}(1) velB(1) velB(2) V{i,j}(4)];
            
            % Calculate south face flux
            s_flux{i,j} = flux(Vb,nB);
        else
            % Check to see is contravariant M >= 1
            if conM >= 1
                % Calculate south face flux using full flux
                s_flux{i,j} = flux(V{i,j},nS{i,j});
            else
                % Calculate south face flux using the negative
                % of the adjacent cell north face flux
                s_flux{i,j} = -1*n_flux{i,j-1};                            
            end
        end
        
        % Assemble residual:
        % SUM(face flux * face area)
        resid{i,j} = e_flux{i,j}*sE(i,j) + ...
                     n_flux{i,j}*sN(i,j) + ...
                     w_flux{i,j}*sW(i,j) + ...
                     s_flux{i,j}*sS(i,j);
    end
end
end

% CALCULATE TIME REMAINING 
function time = fixTime(s)
% Input: time in seconds

min = fix(s/60);
sec = rem(s,60);
sec = fix(sec);

time = sprintf('%d:%02d',min,sec);
end


