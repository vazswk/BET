function [beta,rho] = betas_ternary(X,I,alpha)

% beta ... embedding change probs
% profile ... costs sorted by increasing value
load('ixlnx3.mat');       

payload = alpha*numel(X);     % Absolute payload in nats

L = 10^3;                     % Initial search interval for lambda
R = 10^6;

fL = sum(h_tern(1./invxlnx3_fast(L*I,ixlnx3))) - payload;
fR = sum(h_tern(1./invxlnx3_fast(R*I,ixlnx3))) - payload;

while fL*fR > 0             % If the range [L,R] does not cover alpha,
    if fL > 0, R = 2*R;     % enlarge the search interval
    else       L = L/2; end
    fL = sum(h_tern(1./invxlnx3_fast(L*I,ixlnx3))) - payload;
    fR = sum(h_tern(1./invxlnx3_fast(R*I,ixlnx3))) - payload;
end

i = 0;
fM = 1; 
while ((abs(fM)>0.0001) && (i<60)) %*(L+R)/2--(abs(L-R)>0.00001) && 
    M = (L+R)/2;
    fM = sum(h_tern(1./invxlnx3_fast(M*I,ixlnx3))) - payload;
    if fL*fM < 0, R = M; fR = fM;
    else          L = M; fL = fM; end
    %fprintf('%d    %f     %f    %f\n',i,fM,abs(L-R),M);
    i = i+1;
end

beta = 1./invxlnx3_fast(M*I,ixlnx3);
rho = log(invxlnx3_fast(M*I,ixlnx3) - 2);
%profile = sort(rho)/max(rho);
