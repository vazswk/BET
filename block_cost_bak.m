function [ bc ] = block_cost( entropy )
%BLOCK_COST 此处显示有关此函数的摘要
%   此处显示详细说明
avg_ent = colfilt(entropy, [8 8], 'distinct', @avg_entropy);
changes = zeros(size(avg_ent));
table= generateTable(4000);

[h, w] = size(avg_ent);

for i = 1:8:h
    for j =1:8:w
        [~,a_idx] = min(abs(table(2,:) - avg_ent(i, j)));
        changes(i:i+7, j:j+7) = table(1, a_idx);
    end
end


bc = log(1./changes -2);

    function avg = avg_entropy(entropy)
        sum_entropy = sum(entropy, 1)/64;
        avg = repmat(sum_entropy, 64, 1);
    end


end

