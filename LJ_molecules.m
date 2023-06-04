epsilon = 1;
sigma = 1;

it_max = 1e6;
tol = 2e-3;
tol_star = 7e-6;

dim = 3;
N = 8;

% X corresponds to the vectors x1, ..., x_N representing the positions of
% each atom. Each i^th row of X is the x_i vector
X = sym('x', [N, dim]);

energy = 0;
for i = 1:N
    for j = 1:i
        if i~=j 
            switch dim
              case 1
                   temp=sqrt( (X(i,1)-X(j,1))^2);
              case 2
                   temp=sqrt( (X(i,1)-X(j,1))^2 + (X(i,2)-X(j,2))^2);
              case 3
                  temp=sqrt( (X(i,1)-X(j,1))^2 + (X(i,2)-X(j,2))^2 + (X(i,3)-X(j,3))^2 );
            end
            energy = energy + 4*epsilon*( (sigma/temp)^12 - (sigma/temp)^6 );
        end
    end
end

Diff_E = gradient(energy);
Hess_E = hessian(energy);

x_init = rand(N, dim);

%   For N=4 (tetrahedral)
% x_init = [0.448636593326791,0.506783267055629,0.860180503723948;0.884474714194616,0.0760229085454987,0.634390055837083;0.853718337359078,0.362406884500106,0.731457114777647;0.647736744115042,0.231091155805152,0.156355283373281];

%   For N=4 (square-planar)                                  
% x_init = [0.439352012713857,0.0490852758965976,0.0514125253940096;0.591769136380890,0.719694172744749,0.948004653343517;0.617984610357016,0.740092862384778,0.877386850235409;0.454096672485187,0.404798794595577,0.827649947688969];

%   For  N=5 (square pyramid)
% x_init = [1.5,1.5,0;1.5,0,0;0,1.5,0;0,0,0;0.75,0.75,2.1] ;

%   For  N=5 (triangular bi-pyramid)
% x_init = [0,0,0;1,0,0;0.5,1,0;0.5,0.5,1;0.5,0.5,-1];

%   For  N=6 (square bi-pyramid)
% x_init = [1.5,1.5,0;1.5,0,0;0,1.5,0;0,0,0;0.75,0.75,2.1;0.75,0.75,-2.1] ;

%   For  N=6 (square pyramid w\ barycenter)
% x_init = [0,0,0;2,0,0;0,2,0;2,2,0;1,1,2;1,1,1] ;  

%   For  N=8 (cube)
% x_init = [0,0,0;2,0,0;0,2,0;2,2,0; 0,0,2;2,0,2;0,2,2;2,2,2] ;  
                  

x = solve_NonLinear_CG_it(x_init, N, dim, X, Diff_E, Hess_E, tol, tol_star);
E=double(subs(energy,X,x));
visualize_molecule2(x,3,3,strcat("$E=",num2str(E),"$ ","($N=",num2str(N),"$)"));
set(gca,'fontsize',30)

