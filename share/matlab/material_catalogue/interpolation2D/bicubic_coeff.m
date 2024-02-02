function [c] = bicubic_coeff(E, dEda, dEdb, dEdadb, da, db)
% Berechnung der Interpolationskoeffizienten für ein bestimmtes Intervall
% A=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;...
%   0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;...
%  -3 0 0 3 0 0 0 0 -2 0 0 -1 0 0 0 0;...
% 2 0 0 -2 0 0 0 0 1 0 0 1 0 0 0 0;...
% 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;...
% 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;...
% 0 0 0 0 -3 0 0 3 0 0 0 0 -2 0 0 -1;...
% 0 0 0 0 2 0 0 -2 0 0 0 0 1 0 0 1;...
% -3 3 0 0 -2 -1 0 0 0 0 0 0 0 0 0 0;...
% 0 0 0 0 0 0 0 0 -3 3 0 0 -2 -1 0 0;...
% 9 -9 9 -9 6 3 -3 -6 6 -6 -3 3 4 2 1 2;...
% -6 6 -6 6 -4 -2 2 4 -3 3 3 -3 -2 -1 -1 -2;...
% 2 -2 0 0 1 1 0 0 0 0 0 0 0 0 0 0;...
% 0 0 0 0 0 0 0 0 2 -2 0 0 1 1 0 0;...
% -6 6 -6 6 -3 -3 3 3 -4 4 2 -2 -2 -2 -1 -1;...
% 4 -4 4 -4 2 2 -2 -2 2 -2 -2 2 1 1 1 1];
A = [
     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0;
     0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0;
    -3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0;
     0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0;
     9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1;
    -6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1;
     2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0;
    -6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1;
     4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1
];

dadb=da*db;
x = zeros(16,1);
for i=1:4 
x(i)=E(i);
x(i+4)=dEda(i)*da;
x(i+8)=dEdb(i)*db;
x(i+12)=dEdadb(i)*dadb;
end
c = A*x;
end

