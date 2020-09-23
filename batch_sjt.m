%% mian.m
clc;clear all;
parpool(12)
input = 'E:\suwenkang\data\Q75_mod'; 
output = 'E:\suwenkang\project\SJT\stego\Q75';
QF = 75;
wetConst = 10^13;
feature_path_cover_GFR_and_Integral_ccJRM = 'E:\suwenkang\project\SJT\feature\ccJRM_cover_2000_Q75.mat';
feature_path_cover_GFR   = 'E:\suwenkang\project\SJT\feature\GFR_cover_2000_Q75.mat';
my_GFR_and_Integral_ccJRM(input,feature_path_cover_GFR_and_Integral_ccJRM,QF);
my_GFR(input, feature_path_cover_GFR, QF);

CAPA = [0.3]; 
params = 123; % the secret key 
err_GFR_and_Integral_ccJRM = zeros(1,length(CAPA));
err_GFR = zeros(1,length(CAPA));

for x = 1:length(CAPA)
    rate = CAPA(x);

    Output_path = [output '\' num2str(rate)];
    feature_path_stego_GFR_and_Integral_ccJRM  = ['E:\suwenkang\project\SJT\feature\GFR_and_Integral_ccJRM_stego_2000_stc_Q75_' num2str(rate*100) '.mat'];
    feature_path_stego_GFR    = ['E:\suwenkang\project\SJT\feature\GFR_stego_2000_stc_Q75_' num2str(rate*100) '.mat'];
        
    if exist(Output_path,'dir'); rmdir(Output_path,'s'); end    
    if ~exist(Output_path,'dir'); mkdir(Output_path); end

    flist = dir([input '\*.jpg']);
    flen = length(flist);
    fprintf('%s%d\n', 'the num of the files: ',flen);
      for pic_num = 1:flen/500
         for i = (pic_num-1)*500 +1 :pic_num*500
               fprintf('%s\n', ['      processing stage: ' num2str(pic_num)]);
               fprintf('%d%s\n',i, ['      processing image: ' flist(i).name]);
               in_file_name = [input '\' flist(i).name];
               stego_name = [ Output_path '\' flist(i).name];

               img = jpeg_read(in_file_name);               
               dct_coef = double(img.coef_arrays{1}); 
               dct_coef2 = dct_coef;
               dct_coef2(1:8:end,1:8:end) = 0;
               nz_index = find(dct_coef2 ~=0);
               nz_number = length(nz_index);                                  
               q_tab = img.quant_tables{1}; 
               q_matrix = repmat(q_tab,[64 64]);
               
               spatial_cover = double(imread(in_file_name));
               JPEG_cost = SJT(spatial_cover, q_matrix, nz_number);
               
               
%%%%%%%%%%%%%%%%%%%   embedding   %%%%%%%%%%%%%%%%%%
    %%    STC                              
%                rhoM1 = rho;
%                rhoP1 = rho;
% 
%                rhoM1(rhoM1 > wetConst) = wetConst;
%                rhoM1(isnan(rhoM1)) = wetConst;
%                rhoM1(dct_coef < -1023) = wetConst;
% 
%                rhoP1(rhoP1 > wetConst) = wetConst;
%                rhoP1(isnan(rhoP1)) = wetConst;    
%                rhoP1(dct_coef > 1023) = wetConst;   
% 
%                rand('state',123);  % Pseudo-random Permutation for cover elements
%                r_index = randperm(length(nz_index));
%                nz_dct_coef = dct_coef(nz_index(r_index));
% 
%                rand('state',123); % Pseudo-random Permutation for message
%                hidden_message = double(rand(ceil(rate*nz_number),1)>0.5);
% 
%                costs = zeros(3, length(nz_index), 'single'); % for each pixel, assign cost of being changed
%                costs(1,:) = rhoM1(nz_index(r_index));       % cost of changing the first cover pixel by -1, 0, +1
%                costs(3,:) = rhoP1(nz_index(r_index));       % cost of changing the first cover pixel by -1, 0, +1
% 
%                [d stego n_msg_bits l] = stc_pm1_pls_embed(int32(nz_dct_coef)', costs, uint8(hidden_message)', 10); % ternary STC embedding 
%                em_dct_coef = dct_coef;
%                em_dct_coef(nz_index(r_index)) = stego;
% 
%                S_STRUCT = img;
%                S_STRUCT.coef_arrays{1} = em_dct_coef;  
%                jpeg_write(S_STRUCT, stego_name); 
             %%     simulator       
               H = 0;
        %        [stego, dist] = f_embedding_jpg(dct_coef, JPEG_cost, rate, nz_number, H, params, WET);% Embed the secret message  
               [stego, dist] = f_sim_embedding_jpg_2(dct_coef, JPEG_cost, rate, nz_number, params);
               S_struct = img;
               em_dct_coef = double(stego);  
               S_struct.coef_arrays{1} = em_dct_coef;
               jpeg_write (S_struct,stego_name);               
               
               err_GFR
			   err_GFR_and_Integral_ccJRM
         end
     end     
   %% 
   delete(gcp)
   pause(5)
   
   parpool(12)
   my_GFR(Output_path,feature_path_stego_GFR,QF);
   delete(gcp)
   [test_error, err_std] = my_ensemble_2(feature_path_cover_GFR,feature_path_stego_GFR);
   err_GFR(x) = test_error;
   save('err_GUED_GFR','err_GFR')    
   
   parpool(12)
   my_GFR_and_Integral_ccJRM(Output_path,feature_path_stego_GFR_and_Integral_ccJRM,QF);
   delete(gcp)
   [test_error, err_std] = my_ensemble_2(feature_path_cover_GFR_and_Integral_ccJRM,feature_path_stego_GFR_and_Integral_ccJRM);
   err_GFR_and_Integral_ccJRM(x) = test_error;
   save('err_GUED_GFR_and_Integral_ccJRM','err_GFR_and_Integral_ccJRM')
   
   parpool(12)

end
delete(gcp)


