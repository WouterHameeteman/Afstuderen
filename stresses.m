% Number of fibers
n = 19;

% Define constants for matrix stress
mu = 40e-3; % shear modulus
nu = 0.48; % Poissoin's ratio
kappa = (2*mu*(1+nu))/(3*(1-2*nu)); % compression modulus

% Define constant for collagen fiber stress
k1 = 0.28e-3; % stiffness parameter 1
k2 = 4.75; % stiffness parameter 2
     
% Volume fraction of the matrix
phim = 0.8; 

% Define the stretch ratios, and the time of each fiber
lambdap = [1.2; 1.3; 1.1; 1.1; 1.2; 1.23; 1.24; 1.25; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.2; 1.4];
lambdaptime = randi([1 n],1,n);
lambdae = 1.2;

% Locations of fibers that are still acceptable, that have not overstayed their welcome 
locations = find(lambdaptime > i);

% Volume fraction of collagen fibers, total should equal 1.
phicf = (1 - phim) / length(locations);

%Define some more math symbols
lambda = lambdae*lambdap;
Fe = lambda;
transFe = transpose(Fe);
B = Fe*transFe;
I = 1;
lambda = sqrt(transFe*Fe); % getal vs kolom

% Calculate the J
for m = 1:n
    J(m) = det(Fe(m));
    % stress in the matrix
    sigm(m) = phim*(mu/J(m)*(B(m)-I) + kappa*(log(J(m))/J(m))*I);
    J(m);
end
 
% Create the fibers, and give all of them a number
sigcf = zeros(n,3);
for k = 1:n
    lambdap(k);
    lambdafin = lambdae*lambdap;
    sigcf(k,1) = phicf*(k1*(lambdafin(k))^2)/J(k)*(exp(k2*(lambdafin(k))^2-1)-1);
    sigcf(k,2) = k;
    sigcf(k,3) = lambdafin(k);
end

% Calculate the total stress for each combination of fibers
sumsig = zeros(n);
for i = 1:n
   for j = 1:length(locations)
       sumsig(i) = sumsig(i) + sigcf(j);
   end
   disp('Fibers discontinued:'), find(lambdaptime<=i)
end

% Total stress for matrix + collagen fibers
sigtmatrix = [];
for t = 1:n
    sigt = sigm(t) + sumsig(t);
    sigtmatrix(end+1) = sigt;
   
end
sigtmatrix
scatter(1:n, sigtmatrix);