The VSFs_Nl file is a Matlab/Octave code useful for the computation of the sensitivity factors.


%Analytical computation of the time-independent mean partial pathlengths (MPPLs), mean time sensitivity factors (MTSFs) 
% and variance sensitivity factors (VSFs) of photons in an N-layered semiinfinite optically turbid cylinder. 
% Further details regarding the MPPLs can be found in:

% - http://dx.doi.org/10.1364/BOE.449514
% - http://dx.doi.org/10.7149/OPA.56.2.51145
% - http://dx.doi.org/10.1117/12.2670453

%   Inputs: 
%           ua -> 1-D array of N absorption coefficients [mm^-1]
%           ups -> 1-D array of N reduced scattering coefficients [mm^-1]
%           n -> 1-D array of N refractive indexes [-]
%           l -> 1-D array of N-1 thicknesses [mm]
%           Rad -> Radius of the cylinder [mm]
%           rho -> Array of source-detector distances [mm]
%           Nb -> Number of Bessel zeros [-]
%   Outputs: 
%            MTSFj -> Array of MTSFs [mm s]
%            VSFj -> Array of MTSFs [mm s^2]
%            Lj -> Array of MPPLs [mm]
%            R -> CW reflectance
%            
%
%   Example:
%   N = 10;
%   ua=linspace(0.002,0.01,N); ups=linspace(0.8,1.5,N);
%   n=1.33*ones(1,N); l=2*ones(1,N-1); Rad=500;
%   rho=linspace(0,40,100); Nb = 3000;
%   uaidx=randperm(N); upsidx=randperm(N);
%   ua=ua(uaidx); ups=ups(upsidx);
%   [VSFj,MTSFj,Lj,R] = VSFs_Nl(ua,ups,n,l,Rad,rho,Nb);
%   figure
%   plot(rho,VSFj);
%   figure
%   plot(rho,MTSFj);
%   figure
%   plot(rho,Lj);

% © 2024 by Demián Vera & Héctor García. This work is licensed under CC BY-SA 4.0:
% http://creativecommons.org/licenses/by-sa/4.0/?ref=chooser-v1


-----------------------------------------------------------------------------------

The J0_roots.dat file are the Bessel function of zeroth order's roots. 
The Matlab/Octave file named dK.m is a simple implementation of the Kronecker's delta. 

We recommend to put all these files in the same directory.
