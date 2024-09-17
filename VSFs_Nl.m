function [VSFj, MTSFj, Lj,R] = VSFs_Nl(ua,ups,n,l,Rad,rho,Nb)

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


%Verifications
if length(ua) ~= length(ups) || length(ua) ~= length(n) || length(ua) ~= length(l)+1
    error('The length of the optical parameters arrays must be consistent');
end


Nl = length(ua); % Number of layers

c = 299.7925; % Speed of light in vacuum.
v = c./n; % Speed of light in media.

% Calculation of extrapolated boundary, extrapolated radius and factor A and
% some definitions.

D = 1./(3*ups); %Dispersion coefficient for layer 1 [mm^-1]
n_out = 1;      %Refractive index of external medium (air)
nn = n(1)/n_out;
A_factor = 1.71148-5.27938*nn+5.10161*nn^2-0.66712*nn^3+0.11902*nn^4;
zb = 2*A_factor*D(1);
r_prime = Rad + zb;
z0 = 1/ups(1);


%Reading of Bessel zeros file
bzfile = load('J0_roots.dat','r');
bz = bzfile(1:Nb);
bzeros = bz/r_prime;
sn = bzeros.^2;

% Computation of alpha.
a = sqrt(ua./D + sn);

beta(:,Nl) = D(Nl-1).*a(:,Nl-1).*n(Nl-1).^2.*cosh(a(:,Nl-1).*l(Nl-1)) + D(Nl).*a(:,Nl).*n(Nl).^2.*sinh(a(:,Nl-1).*l(Nl-1)); % Eq. (23).
gamma(:,Nl) = D(Nl-1).*a(:,Nl-1).*n(Nl-1).^2.*sinh(a(:,Nl-1).*l(Nl-1)) + D(Nl).*a(:,Nl).*n(Nl).^2.*cosh(a(:,Nl-1).*l(Nl-1)); % Eq. (24).

if Nl<3 % Two-layered medium
    beta = ones(Nb,Nl+1);
    gamma = ones(Nb,Nl+1);

elseif Nl <= 3
    %3-layered medium

    beta(:,3) = beta(:,end);
    gamma(:,3) = gamma(:,end);

else

    for ii=Nl-1:-1:3
        beta(:,ii) = D(ii-1).*a(:,ii-1).*n(ii-1).^2.*cosh(a(:,ii-1).*l(ii-1)).*beta(:,ii+1)+...
            D(ii).*a(:,ii).*n(ii).^2.*sinh(a(:,ii-1).*l(ii-1)).*gamma(:,ii+1);% Eq. (21).
        gamma(:,ii) = D(ii-1).*a(:,ii-1).*n(ii-1).^2.*sinh(a(:,ii-1).*l(ii-1)).*beta(:,ii+1)+...
            D(ii).*a(:,ii).*n(ii).^2.*cosh(a(:, ii-1).*l(ii-1)).*gamma(:,ii+1);% Eq. (22).
    end

end


% Reflectance
[R,B,C] = reflectance(a,ups,n,l,Rad,rho,bzeros, beta, gamma); % Eq. (2).

% First order derivatives

ap1 = 1./(2.*D.*a); % Eq. (32).


%% Initial values of beta and gamma

fN_1 = D(Nl-1).*n(Nl-1).^2; % Just notation.
fN = D(Nl).*n(Nl).^2;       % Just notation.
arg = a(:,Nl-1)*l(Nl-1); % Argument of the hyperbolic functions of the initial values.

% Allocation.

dB = zeros(Nb,Nl,Nl);
dG = zeros(Nb,Nl,Nl);

dB2 = zeros(Nb,Nl,Nl,Nl);
dG2 = zeros(Nb,Nl,Nl,Nl);

dB3 = zeros(Nb,Nl,Nl,Nl,Nl);
dG3 = zeros(Nb,Nl,Nl,Nl,Nl);

% Third order derivatives (initial values)

for j = Nl-1:Nl

    dB(:,Nl,j) = (fN_1+fN.*a(:,Nl).*l(Nl-1)).*dK(Nl-1,j).*cosh(arg)+...
        (fN_1.*a(:,Nl).*l(Nl-1).*dK(Nl-1,j)+fN.*dK(Nl,j)).*sinh(arg); % Eq. (27).

    dG(:,Nl,j) = (fN_1+fN.*a(:,Nl).*l(Nl-1)).*dK(Nl-1,j).*sinh(arg)+...
        (fN_1.*a(:,Nl).*l(Nl-1).*dK(Nl-1,j)+fN.*dK(Nl,j)).*cosh(arg); % Eq. (28).

    for m = Nl-1:Nl

        dB2(:,Nl,j,m) = cosh(arg).*(fN*l(Nl-1).*(dK(Nl,m).*dK(Nl-1,j)+dK(Nl-1,m).*dK(Nl,j)) +...
            fN_1.*l(Nl-1)^2.*dK(Nl-1,m).*dK(Nl-1,j).*a(:,Nl-1)) + ...
            sinh(arg).*(2*fN_1.*dK(Nl-1,m).*dK(Nl-1,j).*l(Nl-1)+...
            fN.*a(:,Nl).*l(Nl-1)^2.*dK(Nl-1,m).*dK(Nl-1,j)); % Eq. (50).

        dG2(:,Nl,j,m) = sinh(arg).*(fN*l(Nl-1).*(dK(Nl,m).*dK(Nl-1,j)+dK(Nl-1,m).*dK(Nl,j)) +...
            fN_1.*l(Nl-1)^2.*dK(Nl-1,m).*dK(Nl-1,j).*a(:,Nl-1)) + ...
            cosh(arg).*(2*fN_1.*dK(Nl-1,m).*dK(Nl-1,j).*l(Nl-1)+ ...
            fN.*a(:,Nl).*l(Nl-1)^2.*dK(Nl-1,m).*dK(Nl-1,j)); % Eq. (51).

        for nn = Nl-1:Nl

            dB3(:,Nl,j,m,nn) = l(Nl-1)^2.*dK(Nl-1,m).*(cosh(arg).*(fN_1*(2+ dK(Nl-1,j).* dK(Nl-1,nn))+...
                fN.*l(Nl-1).*a(:,Nl)*dK(Nl-1,j) ) + fN_1.*dK(Nl-1,j).*dK(Nl-1,nn) +...
                sinh(arg).*(fN.*(2.*dK(Nl-1,j).*dK(Nl,m)+dK(Nl,j)+fN_1.*l(Nl-1)*a(:,Nl-1).*dK(Nl-1,j)))); % Eq. (67) of the supplementary material.

            dG3(:,Nl,j,m,nn) = l(Nl-1)^2.*dK(Nl-1,m).*(sinh(arg).*(fN_1*(2+ dK(Nl-1,j).* dK(Nl-1,nn))+...
                fN.*l(Nl-1).*a(:,Nl)*dK(Nl-1,j) ) + fN_1.*dK(Nl-1,j).*dK(Nl-1,nn) +...
                cosh(arg).*(fN.*(2.*dK(Nl-1,j).*dK(Nl,m)+dK(Nl,j)+fN_1.*l(Nl-1)*a(:,Nl-1).*dK(Nl-1,j)))); % Eq. (68) of the supplementary material.

        end
    end
end


% Third order derivatives (recurrence relations)

for jj = Nl-1:-1:3
    f1 = D(jj-1).*n(jj-1).^2; % Just notation.
    f2 = D(jj).*n(jj).^2;     % Just notation.
    arg = a(:,jj-1).*l(jj-1); % Just notation.
    for m = 1:Nl

        sqBracket1(:,m) = D(jj-1).*n(jj-1).^2.*(beta(:,jj+1).*dK(jj-1,m) + a(:,jj-1).*dB(:,jj+1,m)) + D(jj).*n(jj).^2.*a(:,jj).*dK(jj-1,m).*gamma(:,jj+1).*l(jj-1) ;

        sqBracket2(:,m) = D(jj).*n(jj).^2.*(dK(jj,m).*gamma(:,jj+1)+dG(:,jj+1,m).*a(:,jj) ) +D(jj-1).*n(jj-1).^2.*dK(jj-1,m).*a(:,jj-1).*beta(:,jj+1).*l(jj-1);

        dB(:,jj,m) = sqBracket1(:,m).*cosh(arg) +sqBracket2(:,m).*sinh(arg) ;

        dG(:,jj,m) =  sqBracket1(:,m).*sinh(arg) +sqBracket2(:,m).*cosh(arg) ;


        for k =1:Nl

            AA = f1.*(dK(jj-1,m).*beta(:,jj+1)+a(:,jj-1).*dB(:,jj+1,m)) + f2.*dK(jj-1,m).*l(jj-1).*a(:,jj).*gamma(:,jj+1);

            BB = f2.*(dK(jj,m).*gamma(:,jj+1)+a(:,jj).*dG(:,jj+1,m)) + f1.*dK(jj-1,m).*l(jj-1).*a(:,jj-1).*beta(:,jj+1);

            AAp = f1.*(dK(jj-1,k).*dB(:,jj+1,m) + dK(jj-1,m).*dB(:,jj+1,k)+a(:,jj-1).*dB2(:,jj+1,m,k)) +...
                f2.*dK(jj-1,m).*l(jj-1).*(dK(jj,k).*gamma(:,jj+1) + a(:,jj).*dG(:,jj+1,k));

            BBp = f2.*(dK(jj,k).*dG(:,jj+1,m) + dK(jj,m).*dG(:,jj+1,k)+a(:,jj).*dG2(:,jj+1,m,k)) +...
                f1.*dK(jj-1,m).*l(jj-1).*(dK(jj-1,k).*beta(:,jj+1) + a(:,jj-1).*dB(:,jj+1,k));

            dB2(:,jj,m,k) = (AAp+BB.*dK(jj-1,k).*l(jj-1)).*cosh(arg) +...
                (BBp+AA.*dK(jj-1,k).*l(jj-1)).*sinh(arg);

            dG2(:,jj,m,k) = (AAp+BB.*dK(jj-1,k).*l(jj-1)).*sinh(arg) +...
                (BBp+AA.*dK(jj-1,k).*l(jj-1)).*cosh(arg);

            dB2(:,jj,k,m) = dB2(:,jj,m,k);
            dG2(:,jj,k,m) = dG2(:,jj,m,k);

            for p = 1:Nl

                AApn = f1.*(dK(jj-1,p).*dB(:,jj+1,m) + dK(jj-1,m).*dB(:,jj+1,p)+a(:,jj-1).*dB2(:,jj+1,m,p)) +...
                    f2.*dK(jj-1,m).*l(jj-1).*(dK(jj,p).*gamma(:,jj+1) + a(:,jj).*dG(:,jj+1,p));

                BBpn = f2.*(dK(jj,p).*dG(:,jj+1,m) + dK(jj,m).*dG(:,jj+1,p)+a(:,jj).*dG2(:,jj+1,m,p)) +...
                    f1.*dK(jj-1,m).*l(jj-1).*(dK(jj-1,p).*beta(:,jj+1) + a(:,jj-1).*dB(:,jj+1,p));

                % Eq. (65).
                dLambda= f1.*(dK(jj-1,k).*dB2(:,jj+1,m,p) + dK(jj-1,m).*dB2(:,jj+1,k,p) + dK(jj-1,p).*dB2(:,jj+1,m,k) + a(:,jj-1).*dB3(:,jj+1,m,k,p)) + ...
                    f2.*dK(jj-1,k).*l(jj-1)*(dK(jj,k).*dG(:,jj+1,p) + dK(jj,p).*dG(:,jj+1,k)+a(:,jj).*dG2(:,jj+1,k,p));

                % Eq. (66).
                dGamma = f2.*(dK(jj,k).*dG2(:,jj+1,m,p) + dK(jj,m).*dG2(:,jj+1,k,p)+dK(jj,p).*dG2(:,jj+1,m,k) + a(:,jj).*dG3(:,jj+1,m,k,p)) + ...
                    f1.*dK(jj-1,k).*l(jj-1)*(dK(jj-1,k).*dB(:,jj+1,p)+dK(jj-1,p).*dB(:,jj+1,k)+a(:,jj-1).*dB2(:,jj+1,k,p));

                % Eq. (63).
                dB3(:,jj,m,k,p) = (dLambda +l(jj-1).*dK(jj-1,k).*BBpn+dK(jj-1,p).*l(jj-1).*(BBp+dK(jj-1,k).*l(jj-1).*AA)).*cosh(arg)+...
                    (dGamma +l(jj-1).*dK(jj-1,k).*AApn+dK(jj-1,p).*l(jj-1).*(AAp+dK(jj-1,k).*l(jj-1).*BB)).*sinh(arg);

                % Eq. (64).
                dG3(:,jj,m,k,p) = (dLambda +l(jj-1).*dK(jj-1,k).*BBpn+dK(jj-1,p).*l(jj-1).*(BBp+dK(jj-1,k).*l(jj-1).*AA)).*sinh(arg)+...
                    (dGamma +l(jj-1).*dK(jj-1,k).*AApn+dK(jj-1,p).*l(jj-1).*(AAp+dK(jj-1,k).*l(jj-1).*BB)).*cosh(arg);

            end
        end
    end
end

% We keep only the value jj = 3.

if Nl==1
    dB =0.*squeeze(dB(:,1,:));
    dG = 0.*squeeze(dG(:,1,:));

    dB2=0.*squeeze(dB2(:,1,:,:));
    dG2=0.*squeeze(dG2(:,1,:,:));

    dB3=0.*squeeze(dB3(:,1,:,:,:));
    dG3=0.*squeeze(dG3(:,1,:,:,:));
elseif Nl ==2
    dB =0.*squeeze(dB(:,1,:));
    dG = 0.*squeeze(dG(:,1,:));

    dB2=0.*squeeze(dB2(:,1,:,:));
    dG2=0.*squeeze(dG2(:,1,:,:));

    dB3=0.*squeeze(dB3(:,1,:,:,:));
    dG3=0.*squeeze(dG3(:,1,:,:,:));
else
    dB =squeeze(dB(:,3,:));
    dG = squeeze(dG(:,3,:));

    dB2=squeeze(dB2(:,3,:,:));
    dG2=squeeze(dG2(:,3,:,:));

    dB3=squeeze(dB3(:,3,:,:,:));
    dG3=squeeze(dG3(:,3,:,:,:));

end

% First layers derivatives (A and B).

S0 = sinh(a(:,1)*(z0+zb));
S1 = sinh(a(:,1)*zb);
S2 =D(:,1).*a(:,1).*exp(a(:,1).*(l(1)+zb));

S0p = (z0+zb).*cosh(a(:,1).*(z0+zb));
S1p = zb.*cosh(a(:,1).*zb);
S2p = D(:,1).*exp(a(:,1).*(l(1)+zb)).*(1+a(:,1).*(l(1)+zb));

Bp = (S1.*S0p + S0.*S1p)./S2 - S0.*S1./S2.^2.*S2p;

x1 = a(:,1).*(z0+2*zb);
x2 = a(:,1).*(z0);

Ap = ((x1+1).*exp(-x1)-(x2+1).*exp(-x2))./(2*D(:,1).*a(:,1).^2);

arg = a(:,1)*(l(1)+zb);

delta = D(1)*n(1)^2.*a(:,1).*beta(:,3) - D(2)*n(2)^2.*a(:,2).*gamma(:,3) ;
Delta = D(1)*n(1)^2.*a(:,1).*beta(:,3).*cosh(arg) + D(2)*n(2)^2.*a(:,2).*gamma(:,3).*sinh(arg) ;

Der = zeros(Nb,Nl);
deltapm = zeros(Nb, Nl);
Deltapm = zeros(Nb, Nl);
Cp = zeros(Nb, Nl);

% Derivative of C

BpCp = zeros(Nb, Nl);

for m =1:Nl
    deltapm(:,m) = D(1).*n(1)^2.*(dK(1,m).*beta(:,3) + dB(:,m).*a(:,1)) - D(2).*n(2)^2.*(dK(2,m).*gamma(:,3)+dG(:,m).*a(:,2)) ;

    Deltapm(:,m) = (D(1)*n(1)^2.*(dK(1,m).*beta(:,3)+dB(:,m).*a(:,1)) + D(2)*n(2)^2*(l(1)+zb).*a(:,2).*gamma(:,3)*dK(1,m)).*cosh(arg) + ...
        (D(1)*n(1)^2*dK(1,m)*(l(1)+zb).*a(:,1).*beta(:,3)+D(2)*n(2)^2.*(dK(2,m).*gamma(:,3)+dG(:,m).*a(:,2))).*sinh(arg);

    Cp(:,m) = (deltapm(:,m).*Delta - Deltapm(:,m).*delta)./Delta.^2;
    BpCp(:,m) = Cp(:,m).*Bp;

    Der(:,m) = Ap.*dK(1,m) + Bp.*C.*dK(1,m) + Cp(:,m).*B;
end

const = (pi*r_prime)^2;

Dert = zeros(Nl, length(rho));
for j = 1:length(rho)
    Dert(:,j) =  squeeze(sum(Der.*ap1.*besselj(0,bzeros.*rho(:,j))./besselj(1,bzeros*r_prime).^2/(const), 1, 'omitnan'));
end

% Mean partial pathlengths
Lj = -Dert./R';


% Second order derivative of A.

x1 = a(:,1).*(z0+2*zb);
x2 = a(:,1).*(z0);


App = -((x1.^2+2.*(x1+1)).*exp(-x1) - (x2.^2+2.*(x2+1)).*exp(-x2))./(2*D(:,1).*a(:,1).^3);
ap2 = -1./(4.*D.^2.*a.^3);

% Second order derivative of B.

S0pp = (z0+zb)^2.*S0;
S1pp = zb^2.*S1;
S2pp = D(1).*(l(1)+zb).*exp(a(:,1).*(l(1)+zb)).*(2+a(:,1).*(l(1)+zb));

term1 = ((S0pp.*S1 + 2.*S1p.*S0p+S1pp.*S0).*S2 - (S0p.*S1+S1p.*S0).*S2p)./S2.^2; % First term of the second order derivative of B.
term2 = (((S0p.*S1+S1p.*S0).*S2p + S0.*S1.*S2pp).*S2.^(-2) - S2.^(-3).*2.*S0.*S1.*S2p.^2);  % Second term of the second order derivative of B.

Bpp = term1 - term2;

% Second order derivative of C

C1 = zeros(Nb, Nl);
C2 = zeros(Nb, Nl);
Cpp = zeros(Nb, Nl,Nl);
deltappkm = zeros(Nb, Nl,Nl);
Deltappkm = zeros(Nb, Nl,Nl);
arg = a(:,1)*(l(1)+zb);

delta = D(1)*n(1)^2.*a(:,1).*beta(:,3) - D(2)*n(2)^2.*a(:,2).*gamma(:,3) ;
Delta = D(1)*n(1)^2.*a(:,1).*beta(:,3).*cosh(arg) + D(2)*n(2)^2.*a(:,2).*gamma(:,3).*sinh(arg) ;

Cp = zeros(Nb, Nl);

% Derivative of C

BpCp = zeros(Nb, Nl);

for m =1:Nl

    Cp(:,m) = (deltapm(:,m).*Delta - Deltapm(:,m).*delta)./Delta.^2;
    BpCp(:,m) = Cp(:,m).*Bp;

end

% Here, the computations needed for Eq. (40) begin, splitted in order to
% ease the calculations. 
for m =1:Nl

    C1(:,m)=D(1)*n(1)^2.*(dK(1,m).*beta(:,3)+a(:,1).*dB(:,m))+D(2)*n(2)^2*a(:,2)*dK(1,m).*(l(1)+zb).*gamma(:,3) ;
    C2(:,m)=D(2)*n(2)^2.*(dK(2,m).*gamma(:,3)+a(:,2).*dG(:,m))+D(1)*n(1)^2*a(:,1)*dK(1,m).*(l(1)+zb).*beta(:,3) ;


    for k =1:Nl

        deltappkm(:,m,k) = D(1).*n(1)^2.*(dK(1,k).*dB(:,m)+dB(:,k).*dK(1,m)+dB2(:,m,k).*a(:,1)) - ...
            D(2)*n(2)^2.*(dK(2,m).*dG(:,k)+dG2(:,m,k) .*a(:,2)+dK(2,k)*dG(:,m));

        D1= (D(1).*n(1)^2.*(dK(1,m).*dB(:,k)+dB(:,m).*dK(1,k)+dB2(:,m,k).*a(:,1))+...
            D(2)*n(2)^2.*(l(1)+zb).*dK(1,m).*(dK(2,m).*gamma(:,3) + dG(:,k).*a(:,2))).*cosh(arg) +...
            dK(1,k).*sinh(arg).*(l(1)+zb).*C1(:,m);

        D2= (D(2).*n(2)^2.*(dK(2,m).*dG(:,k)+dG(:,m).*dK(2,k)+dG2(:,m,k).*a(:,2))+...
            D(1)*n(1)^2.*(l(1)+zb).*dK(1,m).*(dK(1,k).*beta(:,3) + dB(:,k).*a(:,1))).*sinh(arg) +...
            dK(1,k).*cosh(arg).*(l(1)+zb).*C2(:,m);

        Deltappkm(:,m,k)  = D1 + D2;

        Cpp(:,m,k) = deltappkm(:,m,k).*Delta.^(-1) - deltapm(:,m).*Delta.^(-2).*Deltapm(:,k) -...
            Deltappkm(:,m,k).*delta.*Delta.^(-2)-   Deltapm(:,m).*deltapm(:,k).*Delta.^(-2)+...
            2.*Deltapm(:,m).*delta.*Delta.^(-3).*Deltapm(:,k);

        Cpp(:,k,m) = Cpp(:,m,k);

    end
end

d2R_dkm = zeros(Nb,Nl,Nl);
part1 = zeros(Nb,Nl,Nl);
part2 = zeros(Nb,Nl,Nl);
apr2 = repmat(ap2,[1,1,Nl]);
%% Mean time sensitivity factors

for ii = 1:Nl
    for jj =1:Nl

        d2R_dkm(:,ii,jj) = App.*dK(1,ii)*dK(1,jj) + dK(1,ii)*(Bpp.*C.*dK(1,jj) + BpCp(:,jj)) + BpCp(:,ii).*dK(1,jj) +Cpp(:,ii,jj).*B;
        if ii~=jj
            d2R_dkm(:,jj,ii) = d2R_dkm(:,ii,jj);
        end
    end
end

% We apply the chain rule to get the derivative of the Green function with
% respect to the absorption coefficient.
for ii = 1:Nl
    for jj = 1:Nl
        part1(:,ii,jj) = Der(:,ii).*dK(ii,jj).*apr2(:,ii,jj);
        part2(:,ii,jj) = d2R_dkm(:,ii,jj).*ap1(:,jj).*ap1(:,ii);
        if ii~=jj
            part2(:,jj,ii) = part2(:,ii,jj);
        end

    end
end

finalTerm1 = zeros(Nl,Nl,length(rho));
finalTerm2 = zeros(Nl,Nl,length(rho));


BF = besselj(0,bzeros.*rho)./besselj(1,bzeros*r_prime).^2;

for j = 1:length(rho)

    finalTerm1(:,:,j) = squeeze(sum(part1.* BF(:,j)/const, 1, 'omitnan'));
    finalTerm2(:,:,j) = squeeze(sum(part2.* BF(:,j)/const, 1, 'omitnan'));
end

lklm = zeros(Nl,Nl,length(rho));

for i = 1:Nl
    for j =1:Nl
        lklm(i,j,:) = 1./R.*squeeze(finalTerm1(i,j,:) + finalTerm2(i,j,:)); % Cross-term.
        if i ~= j
            lklm(j,i,:) = lklm(i,j,:); % Cross-term.
        end
    end
end

meanTime=squeeze(sum(Lj./v', 1)); % This is why we needed v.

prev1 = zeros(Nl,Nl,length(rho));
prev2 =zeros(Nl,length(rho));

for i = 1:size(lklm,1)
    for j = 1:size(lklm,2)
        prev1(i,j,:) = -squeeze(lklm(i,j,:))./v(j);
        prev2(i,:) =  (meanTime.*Lj(i,:))';
    end
end


% The mean time sensitivity factors.
MTSFj = squeeze(sum(prev1,2)) + prev2;


% Derivatives of C.

% Some notation and allocation, again.

Cppp = zeros(Nb,Nl,Nl,Nl);
L1=(l(1)+zb);
arg=a(:,1).*(l(1)+zb);
f1 = D(1).*n(1).^2;
f2 = D(2).*n(2).^2;

for m = 1:Nl
    for k = 1:Nl
        for p = 1:Nl
            
            % The third order derivative of C has been splitted in various
            % terms in order to make it easy to read.

            delta3 = f1.*(dK(1,k).*dB2(:,m,p)+dK(1,m).*dB2(:,k,p)+dK(1,p).*dB2(:,m,k)+a(:,1).*dB3(:,m,k,p)) -...
                f2.*(dK(2,m).*dG2(:,k,p)+dK(2,k).*dG2(:,m,p)+dK(2,p).*dG2(:,m,k)+a(:,2).*dG3(:,m,k,p));

            C1 = f1.*(dK(1,m).*beta(:,3)+a(:,1).*dB(:,m))+f2.*a(:,2).*gamma(:,3).*L1.*dK(1,m);
            C2 = f2.*(dK(2,m).*gamma(:,3)+a(:,2).*dG(:,m))+f1.*a(:,1).*beta(:,3).*L1.*dK(1,m);

            L12 = f1.*(dK(1,m).*dB(:,k)+dK(1,k).*dB(:,m)+a(:,1).*dB2(:,m,k))+...
                f2.*(dK(2,k).*gamma(:,3)+dG(:,k).*a(:,2)).*L1.*dK(1,m) + C2.*dK(1,k).*L1;

            L22 = f2.*(dK(2,m).*dG(:,k)+dK(2,k).*dG(:,m)+a(:,2).*dG2(:,m,k))+...
                f1.*(dK(1,k).*beta(:,3)+dB(:,k).*a(:,1)).*L1.*dK(1,m) + C1.*dK(1,k).*L1;

            C1p = f1.*(dK(1,m).*dB(:,p)+dB(:,m).*dK(1,p)+dB2(:,m,p).*a(:,1))+...
                f2.*L1.*dK(1,m).*(dK(2,p).*gamma(:,3) + dG(:,p).*a(:,2)) ;

            C2p = f2.*(dK(2,m).*dG(:,p)+dG(:,m).*dK(2,p)+dG2(:,m,p).*a(:,2))+...
                f1.*L1.*dK(1,m).*(dK(1,p).*beta(:,3) + dB(:,p).*a(:,1));

            Delta3_1 = (f1.*(dK(1,m).*dB2(:,k,p) +dK(1,k).*dB2(:,m,p)+dK(1,p).*dB2(:,m,k)+a(:,1).*dB3(:,m,k,p)) + ...
                f2.*L1.*dK(1,m).*(dK(2,k).*dG(:,p)+dK(2,p).*dG(:,k)+a(:,2).*dG2(:,k,p)) + C2p.*dK(1,k).*L1).*cosh(arg)+...
                L12.*dK(1,p).*L1.*sinh(arg);

            Delta3_2 = (f2.*(dK(2,m).*dG2(:,k,p)+ dK(2,k).*dG2(:,m,p)+dK(2,p).*dG2(:,m,k)+a(:,2).*dG3(:,m,k,p))+ ...
                f1.*L1.*dK(1,m).*(dK(1,k).*dB(:,p)+dK(1,p).*dB(:,k)+dB2(:,k,p).*a(:,1))+C1p.*dK(1,k).*L1).*sinh(arg)+...
                L22.*dK(1,p).*L1.*cosh(arg);



            Delta3 = Delta3_1 + Delta3_2;

            term1 = Deltappkm(:,k,p).*deltapm(:,m);
            term2 = deltappkm(:,m,p).*Deltapm(:,k);
            term3 = Deltapm(:,p).*deltappkm(:,m,k);
            term4 = Delta.*delta3;

            term5 = -deltappkm(:,k,p).*Deltapm(:,m);
            term6 = -Deltappkm(:,m,p).*deltapm(:,k);
            term7 = -Delta3.*delta;
            term8 = -Deltappkm(:,m,k).*deltapm(:,p);

            p1k = deltappkm(:,m,k).*Delta + deltapm(:,m).*Deltapm(:,k)-Deltappkm(:,m,k).*delta -deltapm(:,k).*Deltapm(:,m);
            p1p = deltappkm(:,m,p).*Delta + deltapm(:,m).*Deltapm(:,p)-Deltappkm(:,m,p).*delta -deltapm(:,p).*Deltapm(:,m);

            firstPart = (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8).*Delta.^(-2);
            secondPart = 2.*Delta.^(-3).*p1k.*Deltapm(:,p);
            thirdPart= p1p.*(-2.*Delta.^(-3).*Deltapm(:,k));
            fourthPart = (6.*Delta.^(-4).*Deltapm(:,k).*Deltapm(:,p) - 2.*Delta.^(-3).*Deltappkm(:,k,p)).*(deltapm(:,m).*Delta-Deltapm(:,m).*delta);

            
            Cppp(:,m,k,p) = firstPart - secondPart + thirdPart + fourthPart; % Eq. 

            % Symmetrization purposes. (Necessary?).
            Cppp(:,p,m,k) = Cppp(:,m,k,p);
            Cppp(:,k,p,m) = Cppp(:,m,k,p);
            Cppp(:,k,m,p) = Cppp(:,m,k,p);

        end
    end
end

% Some notation.

x1 = a(:,1)*(z0+2*zb);
x2 = a(:,1)*(z0);

Appp = ((x1.^3 + 3*x1.^2 + 6*(x1+1)).*exp(-x1) - (x2.^3 + 3.*x2.^2 + 6.*(x2+1)).*exp(-x2) )./(2*D(1).*a(:,1).^4); % Eq. (54).

Bppp = (exp(-a(:,1).*(l(1) + zb)).*cosh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(z0 + zb).^3)./(D(1).*a(:,1)) - (6.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb))./(D(1).*a(:,1).^4) - (exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(l(1) + zb).^3)./(D(1).*a(:,1)) - (3.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(l(1) + zb).^2)./(D(1).*a(:,1).^2) - (3.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(z0 + zb).^2)./(D(1).*a(:,1).^2) + (zb.^3.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*cosh(a(:,1).*zb))./(D(1).*a(:,1)) - (3.*zb.^2.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb))./(D(1).*a(:,1).^2) + (6.*exp(-a(:,1).*(l(1) + zb)).*cosh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(z0 + zb))./(D(1).*a(:,1).^3) - (6.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(l(1) + zb))./(D(1).*a(:,1).^3) + (6.*zb.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*cosh(a(:,1).*zb))./(D(1).*a(:,1).^3) - (6.*zb.*exp(-a(:,1).*(l(1) + zb)).*cosh(a(:,1).*(z0 + zb)).*cosh(a(:,1).*zb).*(z0 + zb))./(D(1).*a(:,1).^2) + (6.*zb.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*cosh(a(:,1).*zb).*(l(1) + zb))./(D(1).*a(:,1).^2) + (3.*exp(-a(:,1).*(l(1) + zb)).*cosh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(l(1) + zb).^2.*(z0 + zb))./(D(1).*a(:,1)) - (3.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(l(1) + zb).*(z0 + zb).^2)./(D(1).*a(:,1)) + (3.*zb.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*cosh(a(:,1).*zb).*(l(1) + zb).^2)./(D(1).*a(:,1)) + (3.*zb.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*cosh(a(:,1).*zb).*(z0 + zb).^2)./(D(1).*a(:,1)) + (3.*zb.^2.*exp(-a(:,1).*(l(1) + zb)).*cosh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(z0 + zb))./(D(1).*a(:,1)) - (3.*zb.^2.*exp(-a(:,1).*(l(1) + zb)).*sinh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(l(1) + zb))./(D(1).*a(:,1)) + (6.*exp(-a(:,1).*(l(1) + zb)).*cosh(a(:,1).*(z0 + zb)).*sinh(a(:,1).*zb).*(l(1) + zb).*(z0 + zb))./(D(1).*a(:,1).^2) - (6.*zb.*exp(-a(:,1).*(l(1) + zb)).*cosh(a(:,1).*(z0 + zb)).*cosh(a(:,1).*zb).*(l(1) + zb).*(z0 + zb))./(D(1).*a(:,1));

ap3 =3./(8.*D.^3.*a.^5); % Eq. (69).

d3G=zeros(Nb,Nl,Nl,Nl); % Allocation of the third order derivative of the Green function for the first layer.

for m= 1:Nl
    for k=1:Nl
        for p=1:Nl

            % Eq. (53), divided in five parts.
            firstPart = Appp.*dK(1,m).*dK(1,k).*dK(1,p);
            secondPart = Bppp.*dK(1,m).*dK(1,k).*dK(1,p).*C;
            thirdPart = Bpp.*dK(1,k).*dK(1,m).*Cp(:,p) + Bpp.*dK(1,k).*dK(1,p).*Cp(:,m) +  Bpp.*dK(1,m).*dK(1,p).*Cp(:,k);
            fourthPart = Cpp(:,m,k).*Bp.*dK(1,p) + Bp.*dK(1,k).*Cpp(:,m,p) + Bp.*dK(1,m).*Cpp(:,k,p);
            fifthPart = Cppp(:,m,k,p).*B;


            d3G(:,m,k,p) = firstPart + secondPart + thirdPart + fourthPart+ fifthPart;


            d3G(:,k,m,p)=d3G(:,m,k,p);
            d3G(:,p,k,m)=d3G(:,m,k,p);
            d3G(:,m,p,k)=d3G(:,m,k,p);

        end
    end
end

apr3 = repmat(ap3,1,1,Nl,Nl);

firstPart = zeros(Nb,Nl,Nl,Nl);
secondPart = firstPart;
thirdPart = secondPart;

for m = 1:Nl
    for k = 1:Nl
        for p = 1:Nl

            % We apply the chain rule in order to get the derivatives with
            % respect to the absortion coefficient.

            firstPart(:,m,k,p) = d3G(:,m,k,p).*ap1(:,p).*ap1(:,m).*ap1(:,k);

            secondPart(:,m,k,p) = d2R_dkm(:,m,k).*apr2(:,k,p).*dK(k,p).*ap1(:,m) + ...
                d2R_dkm(:,m,k).*apr2(:,m,p).*dK(m,p).*ap1(:,k) +...
                d2R_dkm(:,m,p).*apr2(:,m,k).*dK(k,m).*ap1(:,p);

            thirdPart(:,m,k,p) = Der(:,m).*apr3(:,m,k,p).*dK(m,p).*dK(m,k);


            % Symmetrization purposes.
            firstPart(:,k,m,p)=firstPart(:,m,k,p);
            secondPart(:,k,m,p)=secondPart(:,m,k,p);
            thirdPart(:,k,m,p)=thirdPart(:,m,k,p);

            firstPart(:,p,k,m)=firstPart(:,m,k,p);
            secondPart(:,p,k,m)=secondPart(:,m,k,p);
            thirdPart(:,p,k,m)=thirdPart(:,m,k,p);


            firstPart(:,m,p,k)=firstPart(:,m,k,p);
            secondPart(:,m,p,k)=secondPart(:,m,k,p);
            thirdPart(:,m,p,k)=thirdPart(:,m,k,p);

        end
    end
end

% Transformation of the triple cross-term. 

lmlklp = zeros(Nl,Nl,Nl,length(rho));

for j = 1:length(rho)
    lmlklp(:,:,:,j) = -R(j).^(-1).*squeeze(sum((firstPart+secondPart+thirdPart).* BF(:,j)/const, 1, 'omitnan'));
end


for m = 1:Nl
    for k =1:Nl
        for p=1:Nl
            lmlklp(k,m,p,:)=lmlklp(m,k,p,:);
            lmlklp(p,k,m,:)=lmlklp(m,k,p,:);

            lmlklp(m,p,k,:)=lmlklp(m,k,p,:);
        end
    end
end


meanSquaredTime = squeeze(sum(sum(lklm./(v'.*v),1),2));
meanTime=squeeze(sum(Lj./v', 1));


finalTerm3 = zeros(Nl,Nl,Nl,length(rho));
finalTerm2 = zeros(Nl,Nl,length(rho));
finalTerm1 = zeros(Nl,length(rho));

for i = 1:Nl
    for j = 1:Nl
        for k = 1:Nl

            finalTerm1(i,:) = squeeze(Lj(i,:))'.*(meanSquaredTime - 2.*meanTime'.^2);
            finalTerm2(i,j,:) = 2.*meanTime'.*squeeze(lklm(i,j,:))./v(j);
            finalTerm3(i,j,k,:) = -squeeze(lmlklp(i,j,k,:))./(v(k).*v(j));

        end
    end
end

% Finally, the variance sensitivity factors.
VSFj = squeeze(sum(sum(finalTerm3,2),3)) + squeeze(sum(finalTerm2,2)) + finalTerm1;


end


