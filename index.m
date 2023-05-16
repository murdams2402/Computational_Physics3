% Defining an index function giving the corresponding index of the node's
% currents/voltage when numbered as precised on the problem sheet
function l = index(i, j, N)
    l = N*(i-1) + j;
end