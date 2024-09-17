function [R,B,C]=reflectance(a,ups,n,l,Rad,rho,bzeros, beta,gamma)


D = 1./(3*ups);                      %Dispersion coefficient for layer 1 [mm^-1]
n_out = 1;                          %Refractive index of external medium (air)
nn = n(1)/n_out;
A_factor = 1.71148-5.27938*nn+5.10161*nn^2-0.66712*nn^3+0.11902*nn^4;
zb = 2*A_factor*D(1);
z0 = 1/ups(1);
r_prime = Rad + zb;


exp1 = exp(-a(:,1)*abs(-z0));
exp2 = exp(-a(:,1)*(z0+2*zb));
term1 = (exp1-exp2)./(2*D(1)*a(:,1));


sh1 = sinh(a(:,1)*(z0+zb));
sh2 = sinh(a(:,1)*zb);
exp3 = exp(a(:,1)*(l(1)+zb));
den1 = D(1)*a(:,1).*exp3;
B = sh1.*sh2./den1;


num2 = D(1).*a(:,1).*n(1)^2.*beta(:,3) - D(2).*a(:,2).*n(2)^2.*gamma(:,3);
cacho1 = D(1).*a(:,1).*cosh(a(:,1).*(l(1)+zb)).*n(1)^2.*beta(:,3);
cacho2 = D(2).*a(:,2).*sinh(a(:,1).*(l(1)+zb)).*n(2)^2.*gamma(:,3);
den2 = cacho1 + cacho2;
C = num2./den2;

% Reflectance

green1 = (term1 + B.*C);

const = (pi*r_prime)^2;

BF = besselj(0,bzeros.*rho)./besselj(1,bzeros*r_prime).^2;

R = zeros(length(rho),1);
for j = 1:length(rho)
    R(j) = sum(green1.* BF(:,j))./(const);

end

end