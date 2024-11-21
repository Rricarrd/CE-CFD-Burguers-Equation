function [Ek,k,i,legend,t] = solve_case(N,C1,Ck,Re,LES)

    % Modes and initial velocities
    k = 1:N;
    u_hat_0 = 1./k';
    
    % Solver
    [Ek,i,t] = gauss_seidel(u_hat_0, C1, N, Re, LES, Ck, k);

    % Legend string
    if LES==1
        name = "LES";
    else
        name="DNS";
    end
    
    legend = sprintf("%s with N=%i, Re=%.1f, Ck=%.2g and C1=%.2g",name, N, Re, Ck,C1);


    
end

