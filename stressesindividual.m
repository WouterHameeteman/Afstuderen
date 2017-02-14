function stressesindividual

% [Fe] = calculatelambda;
[n, sigm, sigcf, lambdalijst] = variables();

% totalen optellen om tot totaal van 10 bundels van 10 vezels te komen
t = 1:n;
for i = 1:n
    for j = 1:round(n/2)   
        if sigcf(j,2) == t(i)
            sumsig(j) = sum(sigcf(j:j+9,1));            
        end
    end
end
    
% total stress for matrix + collagen fibers
sigtmatrix = [];
for t = 1:round(n/2)
    sigt = sigm + sumsig(t);
    sigtmatrix(end+1) = sigt
   
end
sigtmatrix;
lambdalijst;
scatter(1:round(n/2), sigtmatrix);
