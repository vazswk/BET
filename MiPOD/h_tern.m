function y = h_tern(x)
% Ternary entropy function expressed in nats
z = x(abs(x-0.5)<0.5-eps); % To prevent underflow
y = -2 * z .* log(z) - (1-2*z) .* log(1-2*z);