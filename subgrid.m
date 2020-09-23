function [ stego_start , pM, p0,pP] = subgrid(coef, se, h_num, w_num, rhoM1, rho0, rhoP1, info_len)
%SUBGRID 此处显示有关此函数的摘要
%   此处显示详细说明
ratio = zeros(h_num, w_num);

[h,w] =size(coef);
sub_len_h = h/h_num;
sub_len_w = w/w_num;
stego_start = zeros(h,w);
pM =  zeros(h,w);
pP =  zeros(h,w);
% infos = zeros(h_num, w_num);
for i=1:h_num
    for j=1:w_num
        h_1 =  1+(i-1)*sub_len_h; h_2 = sub_len_h*i;
        w_1 = 1+(j-1)*sub_len_w;  w_2 = sub_len_w*j;
        
        ratio(i,j) = sum(sum(se( h_1:h_2, w_1:w_2)))/sum(sum(se));
        [stego_start(h_1:h_2, w_1:w_2), lamda, pM(h_1:h_2, w_1:w_2), p0, pP(h_1:h_2, w_1:w_2)]  =  EmbeddingSimulator( coef(h_1:h_2, w_1:w_2), rhoM1(h_1:h_2, w_1:w_2), rho0(h_1:h_2, w_1:w_2), rhoP1(h_1:h_2, w_1:w_2), info_len*ratio(i,j));
    end
end
p0= 1 -pM - pP;

end

