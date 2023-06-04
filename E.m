function [Diff_E, Hess_E] = E(X, N, epsilon, sigma)
    
    function V = U_LJ(r)
        V = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 );
    end

[~, d] = size(X);
energy = 0;

for i = 1:N
    for j = 1:N
        if i~=j 
            switch d
              case 1
                   temp=sqrt( (X(i,1)-X(j,1))^2); %Root to cumpute the distance between atoms
              case 2
                   temp=sqrt( (X(i,1)-X(j,1))^2 + (X(i,2)-X(j,2))^2);
              case 3
                  temp=sqrt( (X(i,1)-X(j,1))^2 + (X(i,2)-X(j,2))^2 + (X(i,3)-X(j,3))^2 );
            end
            energy = energy + U_LJ(temp);
        end
    end
end

Diff_E = gradient(energy);
Hess_E = hessian(energy);
end