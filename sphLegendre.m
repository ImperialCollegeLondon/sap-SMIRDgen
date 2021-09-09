function output = sphLegendre( max_n, X )
% This function calculates the spherical Legendre polynomial function using recursion
% relations over order to reduce computetional complexity
% Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London. based on a technique by Jarret, Habets, and
% Naylor ("Rigid sphere room impulse response simulation: Algorithm and
% applications" 2012 paper)

%   Input(s)
%---------------
%   max_n: maximum spherical harmonic order
%   X:      input of Legendre function

%   output(s)
%---------------
%   output: result of spherical Legendre polynomial function (1 x (max_n+1) -> from 0 to
%   max_nu)
output=zeros(1,max_n+1);
output(1)=1;
output(2)=X;
for n=3:length(output)
    output(n)=(2*n-1)/n * X * output(n-1) - (n-1)/n * output(n-2);
end
end

