function stressesindividual

% [Fe] = calculatelambda;
[n, sigm, sigcf, lambdap, lambdaptime] = variables();

% totalen optellen om tot totaal van 10 bundels van 10 vezels te komen
t = 1:n;
sumsig = zeros(round(n/2));

for i = 1:n
   locations = find(lambdaptime < i);
   for j = 1:length(locations)
       sumsig(i) = sumsig(i) + sigcf(j);
   end
end
        
% total stress for matrix + collagen fibers
sigtmatrix = [];
for t = 1:n
    sigt = sigm(t) + sumsig(t);
    sigtmatrix(end+1) = sigt;
   
end
sigtmatrix
scatter(1:n, sigtmatrix);