function JPEG_cost = SJT(spatial_cover, quant_matrix, num_nzac)
entropy = calc_entropy(spatial_cover, 2*num_nzac);
bc = block_cost(entropy);
JPEG_cost = (bc) .* (quant_matrix);

% entropy = calc_entropy_2(spatial_cover, num_nzac, apha);
% bc = block_cost_2(entropy);
% JPEG_cost = bc.*quant_matrix;

% entropy = calc_entropy_3(spatial_cover, num_nzac, apha, quant_matrix);
% JPEG_cost = block_cost_3(entropy);

end