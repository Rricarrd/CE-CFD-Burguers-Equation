
function [Ek,i,t_elapsed] = gauss_seidel(u_hat_0, C1, N, Re, LES, Ck, k)

% Explicit time integration scheme

% Numerical variables
convergence = 1e-5; % convergence error
error = 1;
i = 0;
dt = C1*Re/N^2; % Time step

% Preallocation
tu = 1; % Initial time unit. Starting from second tu, as tu=1 is known
u_hat(:,tu) = u_hat_0; % Start updating at this point

% Getting time
tic;   

% Iterating until error is 
while error > convergence
    i = i+1;

    % Convective and diffusive terms calculation
	C = convective(u_hat(:,tu),N);
    D = diffusive(u_hat(:,tu),N,Re,Ck,LES);
    dE = D+C; 

    % Next velocity calculation
	u_hat(k,tu+1) = u_hat(k,tu)-dt*dE; % Euler time scheme
	u_hat(1,tu+1) = 1; % Setting u0(t) = 1 to force the solution and instead of Fk
    

    % Error calculation. 
    error = max(abs(u_hat(:,tu+1) - u_hat(:,tu))); % When u is stable finish iterating
    
    % Advance time unit
	tu = tu+1;

end

% Calculating kinetic energy after iterating
Ek(k) = u_hat(k,tu).*conj(u_hat(k,tu));

% Elapsed time
t_elapsed = toc;

end

