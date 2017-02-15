function [n, sigm, sigcf, lambdap] = variables()

% Define constants for matrix stress
mu = 40e-3; % shear modulus
nu = 0.48; % Poissoin's ratio
kappa = (2*mu*(1+nu))/(3*(1-2*nu)); % compression modulus

n = 19;

% Define constant for collagen fiber stress
k1 = 0.28e-3; % stiffness parameter 1
k2 = 4.75; % stiffness parameter 2
     
% Define volume fractions, total should equal 1
phim = 0.8; % matrix
phicf = (1 - phim) / round(n/2); % collagen

lambdap = [1.2; 1.3; 1.1; 1.1; 1.2; 1.23; 1.24; 1.25; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.4];
lambdapr = 1.2;
lambdae = 1.2;
lambda = lambdae*lambdap;
Fe = lambda;

transFe = transpose(Fe);
B = Fe*transFe;
I = 1;
lambda = sqrt(transFe*Fe);

for m = 1:n
    J(m) = det(Fe(m));
    % stress in the matrix
    sigm(m) = phim*(mu/J(m)*(B(m)-I) + kappa*(log(J(m))/J(m))*I);
    J(m);
end
 
%create fibers
sigcf = zeros(n,3);
for k = 1:n
    lambdap(k);
    lambdafin = lambdae*lambdap;
    sigcf(k,1) = phicf*(k1*(lambdafin(k))^2)/J(k)*(exp(k2*(lambdafin(k))^2-1)-1);
    sigcf(k,2) = k;
    sigcf(k,3) = lambdafin(k);
end        