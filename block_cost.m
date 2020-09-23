function [ bc ] = block_cost( entropy )
%BLOCK_COST 此处显示有关此函数的摘要
%   此处显示详细说明
avg_ent = colfilt(entropy, [8 8], 'distinct', @avg_entropy);
table= generateTable(4000);

[h, w] = size(avg_ent);
temp = zeros(int8(h/8), int8(w/8));
for i = 1:8:h
    for j =1:8:w
        [~,a_idx] = min(abs(table(2,:) - avg_ent(i, j)));
        temp((i+7)/8, (j+7)/8) = table(1, a_idx);
    end
end

tempCost = log(1./temp -2);

kernel = fspecial('gaussian',[3 3],1);
tempCost = conv2(tempCost, kernel, 'same');


bc = zeros(h, w);
for i = 1:8:h
    for j =1:8:w
        bc(i:i+7, j:j+7) = tempCost((i+7)/8, (j+7)/8);
    end
end

function avg = avg_entropy(entropy)
    sum_entropy = sum(entropy, 1)/64;
    avg = repmat(sum_entropy, 64, 1);
end


end

