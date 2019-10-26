function [polyArray] = cubspline(x,y)

n = length(x);

%Calculation of interval matrix using array operations
x1 = x; %Create dummy array for substraction to obtain interval widths
x1 = circshift(x1,-1);

h = x1-x; %Calculation of interval matrix
h(n) = []; %Deleting last element of array as doesn't exist in interval matrix

%For the calculations below, D is a (n-2)x(n-2) array
%E is a (n-2)x1 array

h1 = h; %Creating array h1 for use as diagonal in matrix D
h1(n-1) = [];
h1(1) = [];

h2 = circshift(h,-1); %Creating dummy array h2 for calculation of array h3

h3 = 2*(h2 + h); %Array h3 used as diagonal in matrix D
h3(n-1) = [];

A = diag(h1,-1);
B = diag(h1,1);
C = diag(h3,0);

D = A + B + C; %Matrix D created from linear combinated of diagonal matrices

E = zeros(1,n-2);

for j=1:(n-2)
    E(j)=6*(((y(j+2)-y(j+1))/h(j+1))-((y(j+1)-y(j))/h(j))); %Calculation of E matrix values
end

E = E'; %Changing E to a column matrix for calculation

M=inv(D)*E; %This creates a matrix of size (n-2)x1, which contains all m values except m(1) and m(n) which are known (0)

%Adding 2 '0' values and shifting matrix to account for 0 values of M(1) and M(n)
M(n-1) = 0;
M(n) = 0;
M = circshift(M,1);

polyArray = zeros([n-1,4]); %Size initialization of polyArray

for k=1:(n-1)              %Simple calculations using formualas provided
    polyArray(k,4) = y(k);
    polyArray(k,3) = ((y(k+1)-y(k))/h(k))-(h(k)*(2*(M(k))+M(k+1)))/6;
    polyArray(k,2) = M(k)/2;
    polyArray(k,1) = (M(k+1)-M(k))/(6*h(k));
end







