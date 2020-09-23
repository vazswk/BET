function [ cost ] = calc_subcost( coef,q )
%CALC_SUBCOST 此处显示有关此函数的摘要
%   此处显示详细说明
wetcost = 1e13;
block_num = size(coef)/8;
rep_q = repmat(q, block_num(1), block_num(2));

coef = coef./rep_q;
sum_coef = colfilt(coef, [8 8], 'distinct', @sum_blockcoef);

% [h, w] = size(sum_coef);
cost = sum_coef./(abs(coef)+1e-5);
cost(cost>wetcost) = wetcost;

cost(1:8:end, 1:8:end) = (cost(2:8:end, 1:8:end)+cost(1:8:end, 2:8:end))/2;

    function sum_coef = sum_blockcoef(block)
        sum_coef = sum(abs(block), 1);
        sum_coef = repmat(sum_coef, 64, 1);
    end


end

