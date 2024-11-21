
function [D] = diffusive(u_hat,N,Re,Ck,LES)
% Preallocation
D = zeros(N,1);
m = 2;
% Solving for each Fourier mode
for k = 1:N

    if LES == 1 % LES model case
        nu_inf = 0.31*((5-m)/(m+1))*sqrt(3-m)*Ck^(-3/2);
        nu_a = 1 + 34.5*exp(-3.03*N/k); % spectral non_dimensional eddy viscosity
        E_kn = u_hat(N)*conj(u_hat(N));  % Energy cutoff freq

        nu_t = nu_inf*(E_kn/N)^(1/2)*nu_a; % Turbulent viscosity calculated

    else % DNS model case
        nu_t = 0; % No turbulent viscosity
    end
    
    % Diffusive term calculation
    D(k) = k^2*(1/Re + nu_t)*u_hat(k); 
end

end