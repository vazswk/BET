function [y, lambda, pChangeM1, pChange0, pChangeP1] = EmbeddingSimulator_hu_2(x, rhoM1, rhoP1, m)
x = double(x);
n = numel(x);
min_rate = 1e-6;
lambda = calc_lambda_2(rhoP1, rhoM1, m, n);
pChangeP1 = 1/3 - lambda .* rhoP1;
pChangeP1 = max(pChangeP1, min_rate);
pChangeM1 = 1/3 - lambda .* rhoM1;
pChangeM1 = max(pChangeM1, min_rate);
if nargin == 5
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
else
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',sum(100*clock)));
end
randChange = rand(size(x));
y = x;
y(randChange < pChangeP1) = y(randChange < pChangeP1) + 1;
y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) = y(randChange >= pChangeP1 & randChange < pChangeP1+pChangeM1) - 1;
pChange0 = 1 - pChangeP1 - pChangeM1;
end

