function V = U_LJ(r, epsilon, sigma)
if r  > 0
    V = 4*epsilon*( (sigma/r)^12 + (sigma/r)^6 );
else
    V = 0;
end
end