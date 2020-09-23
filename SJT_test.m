%% mian.m
clc;clear all;
currdir = pwd;
input = [currdir '\Cover\'];             % input dir for all QFs
output = [currdir '\Stego\'];        
rate = [0.2];   
params = 123; % the secret key 
WET = single(1e14);% STC wet cost
wetConst = 10^13;

flist = dir([input '\*.jpg']);
flen = length(flist);
fprintf('%s%d\n', 'the num of the files: ',flen);

load('Gaussian_sigmahat_95.mat')  
param_sigmahat = (1./param_sigmahat).^0.8;  % p=[0.7,1.3]

for i=1:length(flist) 
       fprintf('%d%s\n',i, ['      processing image: ' flist(i).name]);
       in_file_name = [input '\' flist(i).name];
       stego_name = [output '\' flist(i).name];
       
       img = jpeg_read(in_file_name);
       dct_coef = double(img.coef_arrays{1}); 
       dct_coef2 = dct_coef; 
       dct_coef2(1:8:end,1:8:end) = 0;
       nz_number = nnz(dct_coef2); % number of non zero ac coefficients         
       q_tab = img.quant_tables{1};
       q_tab(1,1) = 0.5*(q_tab(2,1)+q_tab(1,2));
       q_matrix = repmat(q_tab,[64 64]);
       
%        q_matrix = repmat(param_sigmahat,[64 64]);
       spatial_cover = double(imread(in_file_name));
       JPEG_cost = SJT(spatial_cover, q_matrix, nz_number);
%%%%%%%%%%%%%%%%%%%   embedding   %%%%%%%%%%%%%%%%%%
    %%    STC
%        H = 10;
%        rhoM1 = JPEG_cost;
%        rhoP1 = JPEG_cost;
% 
%        rhoM1(rhoM1 > wetConst) = wetConst;
%        rhoM1(isnan(rhoM1)) = wetConst;
%        rhoM1(dct_coef < -1023) = wetConst;
% 
%        rhoP1(rhoP1 > wetConst) = wetConst;
%        rhoP1(isnan(rhoP1)) = wetConst;    
%        rhoP1(dct_coef > 1023) = wetConst;   
%        
%        rand('state',123);  % Pseudo-random Permutation for cover elements
%        r_index = randperm(length(nz_index));
%        nz_dct_coef = dct_coef(nz_index(r_index));
% 
%        rand('state',123); % Pseudo-random Permutation for message
%        hidden_message = double(rand(ceil(rate*nz_number),1)>0.5);
% 
%        costs = zeros(3, length(nz_index), 'single'); % for each pixel, assign cost of being changed
%        costs(1,:) = rhoM1(nz_index(r_index));       % cost of changing the first cover pixel by -1, 0, +1
%        costs(3,:) = rhoP1(nz_index(r_index));       % cost of changing the first cover pixel by -1, 0, +1
% 
%        [d stego n_msg_bits l] = stc_pm1_pls_embed(int32(nz_dct_coef)', costs, uint8(hidden_message)', H); % ternary STC embedding 
% 
%        em_dct_coef = dct_coef;
%        em_dct_coef(nz_index(r_index)) =stego;
% 
%        S_STRUCT = img;
%        S_STRUCT.coef_arrays{1} = em_dct_coef;

    %%     simulator       
       H = 0;
       [stego, dist] = f_embedding_jpg(dct_coef, JPEG_cost, rate, nz_number, H, params, WET);% Embed the secret message  
%        [stego, dist] = f_sim_embedding_jpg_2(dct_coef, JPEG_cost, rate, nz_number, params);
       S_struct = img;
       em_dct_coef = double(stego);  
       S_struct.coef_arrays{1} = em_dct_coef;
       jpeg_write (S_struct,stego_name);     
end   
D = dct_coef - em_dct_coef;
sum(sum(abs(D(:))))
show_s_dif(dct_coef,em_dct_coef);
J_embed_mod_dif_bar(in_file_name,stego_name)   
[img_h, img_w] = size(D);
for mod_m = 1:8
    for mod_n = 1:8
        D_mode(1+(img_h/8)*(mod_m-1):(img_h/8)*mod_m,1+(img_w/8)*(mod_n-1):(img_w/8)*mod_n) = D(mod_m:8:end,mod_n:8:end);    
    end
end
show_cost_dis(D_mode)
