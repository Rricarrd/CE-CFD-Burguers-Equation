
function [Ck] = convective(u_hat,N)

% Preallocation
Ck = zeros(N,1);

for k=1:N % Solving for each Fourier mode
    
    for q=k-N:N
    p = k-q;

        % If q is negative, get the conjugate of the positive
        if q < 0                            
            u_hat_q = conj(u_hat(-q));
        elseif q > 0                        
            u_hat_q = u_hat(q);
        else
            u_hat_q = 0;
        end  
        
        % If p is negative, get the conjugate of the positive
        if p < 0     
            u_hat_p = conj(u_hat(-p));
        elseif p > 0 
            u_hat_p = u_hat(p);
        else
            u_hat_p = 0;
        end
        
        % Calculation of the convective term
        Ck(k) = Ck(k) + u_hat_p*q*1i*u_hat_q; 


    end   
end

end