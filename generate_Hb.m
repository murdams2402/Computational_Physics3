function [H, b] = generate_Hb(N)
% Inputs:   N: siz of the system
% Outputs:  H: N^2 x N^2 square matrix implementing the poisson equation  
%              with mixed Neumann & Dirichlet boundary conditions 
%           b: solution vector of size N^2 x 1

% Fixing the number of unknown quantities of the system
Nu = N*N;
% N has to be an odd number

% Initialisation
H = zeros(Nu);
b = zeros(Nu, 1);
L=45;
epsilon = 1;
% In this implementation, \Delta x = \Delta y = L/(N-1)
h = L/(N-1);
for i = 1:N
    for j = 1:N
        
        % Implementation of Boundary conditions
        if (i == 1 || i == N) || (j == 1 || j == N)       
            if i == 1 || i == N
                H(index(i, j, N), index(i, j, N)) = 1/h;
                if j >= 1+(N-1)/3 && j <= 1+2*(N-1)/3
                    b(index(i, j, N)) = 1;
                else
                    if i == 1
                        H(index(i, j, N), index(i+1, j, N)) = -1/h;
                    else
                        H(index(i, j, N), index(i-1, j, N)) = -1/h;
                    end
                end

            elseif j == 1 || j == N
                H(index(i, j, N), index(i, j, N)) = 1/h;
                if j == 1
                    H(index(i, j, N), index(i, j+1, N)) = -1/h;
                else
                    H(index(i, j, N), index(i, j-1, N)) = -1/h;
                end
            end
        else
        % Implementation of the Poisson equation
            H(index(i, j, N), index(i, j, N)) = -4/(h*h);
            H(index(i, j, N), index(i-1, j, N)) = 1/(h*h);
            H(index(i, j, N), index(i+1, j, N)) = 1/(h*h);
            H(index(i, j, N), index(i, j-1, N)) = 1/(h*h);
            H(index(i, j, N), index(i, j+1, N)) = 1/(h*h);
            if i == 1+(N-1)/2 && i==j
                b(index(i, j, N)) = -1/epsilon;
            end
        end
    end
end
end