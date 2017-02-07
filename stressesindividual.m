function stressesindividual(Fe)

% Define constants for matrix stress
% This should be changed, because mu and kappa are dependent on eachother
% through Poisson's ratio (nu)
mu = 40e-3; % shear modulus
nu = 0.48; % Poissoin's ratio
kappa = (2*mu*(1+nu))/(3*(1-2*nu)); % compression modulus

% Define constant for collagen fiber stress
k1 = 0.28e-3; % stiffness parameter 1
k2 = 4.75; % stiffness parameter 2
     

        % Define volume fractions, total should equal 1
        phim = 0.8; % matrix
        phicf = 0.02; % collagen
                
        transFe = transpose(Fe);
        J = det(Fe);
        B = Fe*transFe;
        I = 1;
        lambda = sqrt(transFe*Fe);

        % stress in the matrix
        sigm = phim*(mu/J*(B-I) + kappa*(log(J)/J)*I);
        
        %create fibers
        n = 19;
        sigcf = zeros(2,n);
        for k = 1:n
            sigcf(k,1) = phicf*(k1*(lambda)^2)/J*(exp(k2*(lambda)^2-1)-1);
            sigcf(k,2) = k;
        end        
                   
% totalen optellen om tot totaal van 10 bundels van 10 vezels te komen
sumsig = zeros(round(n/2));
% for l = 1:round(n/2)
%     for m = 10 : n+10
%         if sigcf(l,2) < m 
%             sumsig(l) = sum(sigcf(l:l+9,1)); %hij pakt nu nog steeds
%             eerste 10 fibers, en niet fibers die voldoen aan de eis dat
%             ze jonger zijn dan tijd t
%         end
%     end
% end

for l = 1:round(n/2)
    for m = 10 : n+10
        if sigcf(l,2) < m 
            sumsig(l) = 
        end
    end
end
    
% total stress for matrix + collagen fibers
sigtmatrix = [];
for t = 1:round(n/2)
    sigt = sigm + sumsig(t)
    sigtmatrix(end+1) = sigt;
end

scatter(1:round(n/2), sigtmatrix);
