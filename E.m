function [Diff_E, Hess_E] = E(X, N, epsilon, sigma)
    
    function V = U_LJ(r)
        V = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 );
    end

[d, ~] = size(X);
energy = 0;

for i = 1:N
    for j = 1:N
        if i~=j 
            energy = energy + U_LJ(norm(X(1:d, i) - X(1:d, j)));
        end
    end
end

Diff_E = gradient(energy);
Hess_E = hessian(energy);
end