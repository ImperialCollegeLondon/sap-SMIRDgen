function output = sphBesselh( max_nu, Z )
% This function calculates the spherical Hankel function using recursion
% relations over order to reduce computetional complexity
% Implemented by Sina Hafezi. Nov 2014. EEE Dept. Imperial College London. based on a technique by Jarret, Habets, and
% Naylor ("Rigid sphere room impulse response simulation: Algorithm and
% applications" 2012 paper)

%   Input(s)
%---------------
%   max_nu: maximum spherical harmonic order
%   Z:      input of Bessel function

%   output(s)
%---------------
%   output: result of spherical bessel function (1 x (max_nu+1) -> from 0 to
%   max_nu)


output=zeros(1,max_nu+1);
output(1)=exp(i*Z)/(i*Z);
output(2)= -i * (-i/Z + 1/(Z*Z)) * exp(i *Z);
for nu=3:length(output)
    output(nu)=(2*nu-1)/Z * output(nu-1) - output(nu-2);
end

end

