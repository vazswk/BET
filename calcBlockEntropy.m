function [ be ] = calcBlockEntropy( entropy )
%CALCBLOCKENTROPY �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
sum_ent = colfilt(entropy, [8 8], 'distinct', @sum_entropy);

be = sum_ent(1:8:end, 1:8:end);
function avg = sum_entropy(entropy)
        sum_entropy = sum(entropy, 1);
        avg = repmat(sum_entropy, 64, 1);
    end

end

